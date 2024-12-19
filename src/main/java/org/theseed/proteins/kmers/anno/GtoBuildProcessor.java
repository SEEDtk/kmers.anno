/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseMultiReportProcessor;
import org.theseed.genome.Feature;

/**
 * This command is used to apply new annotations and protein families to old genomes in the BV-BRC. We will use
 * a genome source to load the desired genomes, erase the protein families and gene names, and
 * then update the annotations, gene names, and local protein families from the input files. The input files are
 * in a single directory and consist of the following files, all of which are tab-delimited without headers.
 *
 * calls								new annotations (column 2) for each feature ID (column 1)
 * local.family.members.expanded		protein family index number (column 1) and gene name (column 5) for each feature ID (column 2)
 * local.family.defs					function (column 2) for each protein family index number (column 1)
 *
 * An ID will be computed for each family, consisting of the prefix "PLF_", the genus ID, another
 * underscore ("_"), and an eight-digit index number.
 *
 * The positional parameters are the genus number to use in the pseudo-IDs, the name of the input directory
 * containing the annotation files, and the name of the input genome source (file or directory)
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -D	output directory name (default is "gtos" in the current directory)
 * -t	type of genome source (default PATRIC)
 *
 * --clear		erase the output directory before processing
 *
 * @author Bruce Parrello
 *
 */
public class GtoBuildProcessor extends BaseMultiReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GtoBuildProcessor.class);
    /** map of genome IDs to GTOs */
    private Map<String, Genome> genomeMap;
    /** map of protein family IDs to functions */
    private Map<String, String> familyMap;
    /** protein family ID prefix */
    private String genusIdPrefix;
    /** name of annotation file */
    private File annoFile;
    /** name of protein family file */
    private File familyFile;
    /** name of family function file */
    private File functionFile;
    /** input genome source */
    private GenomeSource genomes;
    /** genus ID validation pattern */
    private static final Pattern GENUS_ID_PATTERN = Pattern.compile("[1-9][0-9]*");

    // COMMAND-LINE OPTIONS

    /** type of genome source */
    @Option(name = "--type", aliases = { "--source", "-t" }, usage = "type of input genome source")
    private GenomeSource.Type sourceType;

    /** genus ID for the input genomes */
    @Argument(index = 0, metaVar = "genus_id", usage = "numeric genus ID for the input genomes", required = true)
    private String genusId;

    /** input directory name */
    @Argument(index = 1, metaVar = "inDir", usage = "input directory name for protein family / annotation files", required = true)
    private File inDir;

    /** input genome source (file or directory) */
    @Argument(index = 2, metaVar = "genomeDir", usage = "input genome source (file or directory")
    private File genomeDir;

    @Override
    protected File setDefaultOutputDir(File curDir) {
        return new File(curDir, "gtos");
    }

    @Override
    protected void setMultiReportDefaults() {
        this.sourceType = GenomeSource.Type.PATRIC;
    }

    @Override
    protected void validateMultiReportParms() throws IOException, ParseFailureException {
        // Verify that the genus ID is numeric.
        if (! GENUS_ID_PATTERN.matcher(this.genusId).matches())
            throw new ParseFailureException("Genus ID of \"" + this.genusId + "\" is not valid.");
        // Form the family ID prefix.
        this.genusIdPrefix = "PLF_" + this.genusId + "_";
        // Verify that the input directory exists.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        // Verify that we can read all the input files.
        this.annoFile = new File(this.inDir, "calls");
        if (! this.annoFile.canRead())
            throw new FileNotFoundException("Cannot find the annotation file " + this.annoFile + " in the input directory.");
        this.familyFile = new File(this.inDir, "local.family.members.expanded");
        if (! this.familyFile.canRead())
            throw new FileNotFoundException("Cannot find the family list file " + this.familyFile + " in the input directory.");
        this.functionFile = new File(this.inDir, "local.family.defs");
        if (! this.functionFile.canRead())
            throw new FileNotFoundException("Cannot find the family definition file " + this.functionFile + " in the input directory.");
        // Connect to the genome source.
        this.genomes = this.sourceType.create(this.genomeDir);
        log.info("{} genomes found in source {}.", this.genomes.size(), this.genomeDir);
        // Create the genome map and family map.
        this.genomeMap = new HashMap<String, Genome>(this.genomes.size() * 3 / 2 + 1);
        this.familyMap = new HashMap<String, String>();
    }

    @Override
    protected void runMultiReports() throws Exception {
        // Loop through the genome source, loading it into memory. For each genome, we clear the gene
        // name, annotation, and protein family data for each feature.
        int gCount = 0;
        int pCount = 0;
        for (Genome genome : this.genomes) {
            gCount++;
            log.info("Processing genome #{}: {}", gCount, genome);
            String genomeId = genome.getId();
            // Clean the genome. We need to erase the gene names, the protein families, and the annotated functions.
            for (Feature feat : genome.getFeatures()) {
                if (feat.isProtein()) {
                    feat.setFunction("hypothetical protein");
                    feat.setPgfam(null);
                    feat.setPlfam(null);
                    feat.setGeneName("");
                    pCount++;
                }
            }
            // Put the genome in the map.
            this.genomeMap.put(genomeId, genome);
        }
        log.info("{} genomes read, {} proteins cleared.", gCount, pCount);
        // Read the annotation file and update the annotations.
        try (TabbedLineReader annoStream = new TabbedLineReader(this.annoFile, 4)) {
            // Count the features annotated and the number of bad feature IDs.
            int aCount = 0;
            int errCount = 0;
            log.info("Applying annotations from {}.", this.annoFile);
            long lastMsg = System.currentTimeMillis();
            // Loop through the annotations.
            for (var line : annoStream) {
                // Get the feature to annotate.
                String fid = line.get(0);
                Feature feat = this.getFeature(fid);
                // Insure the feature ID is valid.
                if (feat == null)
                    errCount++;
                else {
                    String function = line.get(1);
                    feat.setFunction(function);
                    aCount++;
                }
                long now = System.currentTimeMillis();
                if (now - lastMsg >= 10000) {
                    log.info("{} features annotated, {} errors.", aCount, errCount);
                    lastMsg = now;
                }
            }
            log.info("{} total features annotated, {} total errors.", aCount, errCount);
        }
        // Read the family definitions. These are used to annotate any proteins that are in families.
        try (TabbedLineReader defStream = new TabbedLineReader(this.functionFile, 6)) {
            log.info("Saving protein family functions from {}.", this.functionFile);
            for (var line : defStream) {
                // Get the family index and function.
                String famIdx = line.get(0);
                String function = line.get(1);
                String famId = this.familyId(famIdx);
                this.familyMap.put(famId,  function);
            }
            log.info("{} family definitions read from {}.", this.familyMap.size(), this.functionFile);
        }
        // Read the family file and update the gene names and protein family definitions.
        try (TabbedLineReader famStream = new TabbedLineReader(this.familyFile, 5)) {
            // Count the families stored, the gene names stored, and the number of bad feature IDs.
            int fCount = 0;
            gCount = 0;
            int errCount = 0;
            int funCount = 0;
            log.info("Applying protein families from {}.", this.familyFile);
            long lastMsg = System.currentTimeMillis();
            // Loop through the family specifications.
            for (var line : famStream) {
                // Get the feature ID and the identified feature.
                String fid = line.get(1);
                Feature feat = this.getFeature(fid);
                if (feat == null)
                    errCount++;
                else {
                    // Form the protein family ID.
                    String famIdx = line.get(0);
                    String plfam = this.familyId(famIdx);
                    // Store the local family ID.
                    feat.setPlfam(plfam);
                    fCount++;
                    // Store the family function.
                    String function = this.familyMap.get(plfam);
                    if (function != null) {
                        feat.setFunction(function);
                        funCount++;
                    }
                    // Check for a gene name.
                    String gene = line.get(4);
                    if (! StringUtils.isBlank(gene)) {
                        feat.setGeneName(gene);
                        gCount++;
                    }
                }
                long now = System.currentTimeMillis();
                if (now - lastMsg >= 10000) {
                    log.info("{} families updated, {} gene names stored, {} functions stored, {} errors.", fCount, gCount, funCount, errCount);
                    lastMsg = now;
                }
            }
            log.info("{} total families updated, {} total gene names stored, {} total functions stored, {} total errors.", fCount, gCount, funCount, errCount);
        }
        // Finally, we save the GTOs to the output directory.
        for (Genome genome : this.genomeMap.values()) {
            File gFile = this.getOutFile(genome.getId() + ".gto");
            log.info("Saving {} to {}.", genome, gFile);
            genome.save(gFile);
        }
    }

    /**
     * @return the full protein family ID for a family index
     *
     * @param famIdx	family index
     */
    private String familyId(String famIdx) {
        return this.genusIdPrefix + StringUtils.leftPad(famIdx, 8, '0');
    }

    /**
     * Locate the specified feature in the GTO map.
     *
     * @param fid	ID of the feature to locate
     *
     * @return the specified feature, or NULL if it is not found
     */
    private Feature getFeature(String fid) {
        Feature retVal = null;
        // Get the genome ID.
        String genomeId = Feature.genomeOf(fid);
        Genome genome = this.genomeMap.get(genomeId);
        // If the genome exists, get the feature itself.
        if (genome != null)
            retVal = genome.getFeature(fid);
        return retVal;
    }

}
