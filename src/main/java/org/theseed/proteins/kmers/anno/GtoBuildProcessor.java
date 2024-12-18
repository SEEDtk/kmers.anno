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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Genome;
import org.theseed.utils.BaseMultiReportProcessor;
import org.theseed.genome.Feature;

/**
 * This command is used to apply new annotations and protein families to old genomes in the BV-BRC. We will use
 * the standard P3Genome object to load genomes from the SOLR data, erase the protein families and gene names, and
 * then update the annotations, gene names, and local protein families from the input files. The input files are
 * in a single directory and consist of the following files, all of which are tab-delimited without headers.
 *
 * calls					new annotations (column 2) for each feature ID (column 1)
 * genome.names				list of genomes to process (column 1)
 * local.family.members		protein family index number (column 1) and gene name (column 5) for each feature ID (column 2)
 *
 * An ID will be computed for each family, consisting of the prefix "PLF_", the genus ID, another
 * underscore ("_"), and an eight-digit index number.
 *
 * The positional parameters are the genus number to use in the pseudo-IDs and the name of the input directory
 * containing the annotation files.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -D	output directory name (default is "gtos" in the current directory)
 *
 * --clear	erase the output directory before processing
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
    /** protein family ID prefix */
    private String genusIdPrefix;
    /** name of annotation file */
    private File annoFile;
    /** name of genome name file */
    private File genomeFile;
    /** name of protein family file */
    private File familyFile;
    /** connection to PATRIC */
    private P3Connection p3;
    /** genus ID validation pattern */
    private static final Pattern GENUS_ID_PATTERN = Pattern.compile("[1-9][0-9]*");

    // COMMAND-LINE OPTIONS

    /** genus ID for the input genomes */
    @Argument(index = 0, metaVar = "genus_id", usage = "numeric genus ID for the input genomes", required = true)
    private String genusId;

    /** input directory name */
    @Argument(index = 1, metaVar = "inDir", usage = "input directory name for protein family / annotation files", required = true)
    private File inDir;

    @Override
    protected File setDefaultOutputDir(File curDir) {
        return new File(curDir, "gtos");
    }

    @Override
    protected void setMultiReportDefaults() {
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
        this.genomeFile = new File(this.inDir, "genome.names");
        if (! this.genomeFile.canRead())
            throw new FileNotFoundException("Cannot find the genome name file " + this.genomeFile + " in the input directory.");
        this.familyFile = new File(this.inDir, "local.family.members");
        if (! this.familyFile.canRead())
            throw new FileNotFoundException("Cannot find the family list file " + this.familyFile + " in the input directory.");
        // Connect to PATRIC.
        log.info("Connecting to PATRIC.");
        this.p3 = new P3Connection();
        // Create the genome map.
        this.genomeMap = new HashMap<String, Genome>();
    }

    @Override
    protected void runMultiReports() throws Exception {
        // Loop through the genome file, reading the GTOs. For each GTO, we clear the gene name, annotation, and protein family data
        // for each feature.
        try (TabbedLineReader genomeStream = new TabbedLineReader(this.genomeFile, 2)) {
            log.info("Reading genomes from {}.", this.genomeFile);
            // Count the genomes read and the proteins erased.
            int gCount = 0;
            int pCount = 0;
            for (var line : genomeStream) {
                String genomeId = line.get(0);
                String genomeName = line.get(1);
                gCount++;
                log.info("Reading genome #{} {}: {}", gCount, genomeId, genomeName);
                Genome genome = P3Genome.load(this.p3, genomeId, P3Genome.Details.FULL);
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
        }
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
        // Read the family file and update the gene names and protein family definitions.
        try (TabbedLineReader famStream = new TabbedLineReader(this.familyFile, 5)) {
            // Count the families stored, the gene names stored, and the number of bad feature IDs.
            int fCount = 0;
            int gCount = 0;
            int errCount = 0;
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
                    String plfam = this.genusIdPrefix + StringUtils.leftPad(famIdx, 8, '0');
                    feat.setPlfam(plfam);
                    fCount++;
                    // Check for a gene name.
                    String gene = line.get(4);
                    if (! StringUtils.isBlank(gene)) {
                        feat.setGeneName(gene);
                        gCount++;
                    }
                }
                long now = System.currentTimeMillis();
                if (now - lastMsg >= 10000) {
                    log.info("{} families updated, {} gene names stored, {} errors.", fCount, gCount, errCount);
                    lastMsg = now;
                }
            }
            log.info("{} total families updated, {} total gene names stored, {} total errors.", fCount, gCount, errCount);
        }
        // Finally, we save the GTOs to the output directory.
        for (Genome genome : this.genomeMap.values()) {
            File gFile = this.getOutFile(genome.getId() + ".gto");
            log.info("Saving {} to {}.", genome, gFile);
            genome.save(gFile);
        }
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
