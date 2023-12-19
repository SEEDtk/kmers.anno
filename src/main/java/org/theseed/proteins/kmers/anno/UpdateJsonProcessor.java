package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;
import java.util.regex.Pattern;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.JsonListInputStream;

/**
 * This sub-command will update the BV-BRC JSON dump files with new annotations in GTOs for the
 * same genomes.  This is a highly specialized command.  For each genome sub-directory with the
 * JSON dump data in it, some files will be copied verbatim, some will require minor field
 * substitutions, and some will need to be completely rebuilt.  Our basic strategy will be to
 * copy the verbatim files first.  Then we will loop through genome_feature.json, replacing the
 * "product" field for each feature.  In parallel, we will rebuild an entry in the subsystem.json
 * file for each feature with subsystems.
 *
 * The positional parameters are the name of the master directory for the JSON genome dumps,
 * the name of the input genome source, and the name of the master directory for the output
 * JSON genome dumps.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -t	type of input genome source (default DIR)
 *
 * --clear	erase the output directory before processing
 */
public class UpdateJsonProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(UpdateJsonProcessor.class);
    /** input genome source */
    private GenomeSource genomes;
    /** input JSON directories */
    private File[] genomeDirs;
    /** set of file names for straight copy */
    private static final Set<String> COPY_FILES = Set.of("genome.json", "protein_structure.json", "sp_gene.json");
    /** genome ID string pattern */
    private static final Pattern GENOME_PATTERN = Pattern.compile("\\d+\\.\\d+");
    /** file filter for subdirectories of the input JSON dump that are genome IDs */
    protected static final FileFilter JGENOME_DIR_FILTER = new FileFilter() {

        @Override
        public boolean accept(File pathname) {
            // Insure we are a directory.
            boolean retVal = pathname.isDirectory();
            if (retVal) {
                String dirName = pathname.getName();
                retVal = GENOME_PATTERN.matcher(dirName).matches();
            }
            return retVal;
        }

    };

    // COMMAND-LINE OPTIONS

    /** if specified, the output directory will be erased before processingg */
    @Option(name = "--clear", usage = "erase the output directory before processing")
    private boolean clearFlag;

    /** input genome source type */
    @Option(name = "--type", aliases = {"-t" }, usage = "input genome source type")
    private GenomeSource.Type sourceType;

    /** input JSON master directory */
    @Argument(index = 0, metaVar = "jsonInDir", usage = "JSON dump input master directory")
    private File jsonInDir;

    /** input genome source */
    @Argument(index = 1, metaVar = "genomeInDir", usage = "input genome source with new annotation data")
    private File genomeInDir;

    /** output directort */
    @Argument(index = 2, metaVar = "jsonOutDir", usage = "JSON dump output master directory")
    private File jsonOutDir;

    @Override
    protected void setDefaults() {
        this.clearFlag = false;
        this.sourceType = GenomeSource.Type.DIR;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Insure the input directory is valid.
        if (! this.jsonInDir.isDirectory())
            throw new FileNotFoundException("Input JSON directory " + this.jsonInDir
                    + " is not found or invalid.");
        log.info("Scanning input directory {}.", this.jsonInDir);
        this.genomeDirs = this.jsonInDir.listFiles(JGENOME_DIR_FILTER);
        if (this.genomeDirs.length <= 0)
            throw new IOException("No genome subdirectories found in " + this.jsonInDir + ".");
        // Get the genome source and verify the genomes in it.
        log.info("Connecting to {} genome source in {},", this.sourceType, this.genomeInDir);
        this.genomes = this.sourceType.create(this.genomeInDir);
        // Insure all the input genomes are in the source.
        int badGenomes = 0;
        Set<String> genomeIds = this.genomes.getIDs();
        for (File gDir : this.genomeDirs) {
            String genomeId = gDir.getName();
            if (! genomeIds.contains(genomeId)) {
                log.error("Missing genome {} in {}.", genomeId, this.genomeInDir);
                badGenomes++;
            }
            throw new ParseFailureException(Integer.toString(badGenomes) + " genomes from "
                    + this.jsonInDir + " not found in " + this.genomeInDir + ".");
        }
        // Set up the output directory.
        if (! this.jsonOutDir.isDirectory()) {
            log.info("Creating output directory {}.", this.jsonOutDir);
            FileUtils.forceMkdir(this.jsonOutDir);
        } else if (this.clearFlag) {
            log.info("Erasing output directory {}.", this.jsonOutDir);
            FileUtils.cleanDirectory(this.jsonOutDir);
        } else
            log.info("Output wll be to {},", this.jsonOutDir);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        final int nGenomes = this.genomeDirs.length;
        log.info("{} genomes to process.", nGenomes);
        // Set up some counters.
        int gCount = 0;
        int substitutions = 0;
        int subRecords = 0;
        int copies = 0;
        // Loop through the input genome JSON dump directories.
        for (File genomeDir : this.genomeDirs) {
            // Get the GTO.
            String genomeId = genomeDir.getName();
            Genome genome = this.genomes.getGenome(genomeId);
            gCount++;
            log.info("Processing genome {} of {}: {}.", gCount, nGenomes, genome);
            // Create the output directory.
            File outDir = new File(this.jsonOutDir, genomeId);
            if (! outDir.isDirectory()) {
                log.info("Creating output directory {}.", outDir);
                FileUtils.forceMkdir(outDir);
            }
            // Copy all the verbatim files.
            for (String name : COPY_FILES) {
                File inFile = new File(genomeDir, name);
                File outFile = new File(outDir, name);
                if (inFile.exists()) {
                    FileUtils.copyFile(inFile, outFile);
                    copies++;
                }
            }
            // Get the input feature file name, and the output feature and subsystem
            // file names.
            File inFile =  new File(genomeDir, "genome_feature.json");
            File outFile = new File(outDir, "genome_feature.json");
            File subFile = new File(outDir, "subsystem.json");
            try (	JsonListInputStream inStream = new JsonListInputStream(inFile);
                    PrintWriter outStream = new PrintWriter(outFile);
                    PrintWriter subStream = new PrintWriter(subFile)) {
                // TODO process copy of the genome
            }
        }
        log.info("{} genomes processed, {} files copied, {} substitutions, {} subsystem records output.",
                gCount, copies, substitutions, subRecords);
    }

}
