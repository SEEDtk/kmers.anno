/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.genome.iterator.GenomeTargetType;
import org.theseed.genome.iterator.IGenomeTarget;
import org.theseed.io.TabbedLineReader;

/**
 * This sub-command applies annotation files to input genomes and writes the results to a new
 * directory.  The annotation files are expected to be those produced by the HashAnnotationProcessor.
 * These are kept in a directory, one file per input genome, with the name XXXXXX.X.anno.tbl,
 * each containing a feature ID (fid), a score (score), and a new annotation (new_annotation).
 *
 * The genomes are input from a genome source and output to a genome target.  Note that some
 * types of targets do not include annotations (LIST, DNAFASTA), so these make no sense.
 *
 * The positional parameters are the name of the annotations directory, the name of the genome
 * source (file or directory), and the name of the genome target (file or directory).  The
 * command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * --source		type of genome source (default DIR)
 * --target		type of genome target (default DIR)
 * --clear		erase the genome target before processing
 */
public class ApplyAnnotationProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ApplyAnnotationProcessor.class);
    /** genome source */
    private GenomeSource genomesIn;
    /** genome target */
    private IGenomeTarget genomesOut;
    /** map of genome IDs to annotation files */
    private Map<String, File> annoMap;
    /** annotation file name pattern */
    private static final Pattern ANNO_FILE_NAME = Pattern.compile("(\\d+\\.\\d+)\\.anno\\.tbl");

    // COMMAND-LINE OPTONS

    /** type of genome source */
    @Option(name = "--source", usage = "type of genome input source")
    private GenomeSource.Type sourceType;

    /** type of genome target */
    @Option(name = "--target", usage = "type of genome output target")
    private GenomeTargetType targetType;

    /** erase-target flag */
    @Option(name = "--clear", usage = "if specified, the genome target will be erased before processing")
    private boolean clearFlag;

    /** annotation directory */
    @Argument(index = 0, metaVar = "annoDir", usage = "name of the annotation file directory", required = true)
    private File annoDir;

    /** input genome source */
    @Argument(index = 1, metaVar = "inDir", usage = "genome source input directory or file", required = true)
    private File inDir;

    /** output genome target */
    @Argument(index = 2, metaVar = "outDir", usage = "genome target output directory or file", required = true)
    private File outDir;

    @Override
    protected void setDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
        this.targetType = GenomeTargetType.DIR;
        this.clearFlag = false;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify the annotations exist.
        if (! this.annoDir.isDirectory())
            throw new FileNotFoundException("Annotation directory " + this.annoDir + " is not found or invalid.");
        File[] annoFiles = this.annoDir.listFiles();
        // Build a map of the annotation files.
        this.annoMap = new TreeMap<String, File>();
        for (File file : annoFiles) {
            if (file.isFile()) {
                // Use the pattern to determine if this is an annotation file.
                String baseName = file.getName();
                Matcher m = ANNO_FILE_NAME.matcher(baseName);
                if (m.matches()) {
                    // Here we have one.
                    if (! file.canRead())
                        throw new IOException("Annotation file " + file + " is unreadable.");
                    // Match group 1 is the genome ID.
                    this.annoMap.put(m.group(1), file);
                }
            }
        }
        log.info("{} annotation files found in {}.", this.annoMap.size(), this.annoDir);
        // Now connect to the genome source.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Input genome source " + this.inDir + " does not exist.");
        log.info("Connecting to {} genome source {}.", this.sourceType, this.inDir);
        this.genomesIn = this.sourceType.create(this.inDir);
        log.info("{} input genomes found.", this.genomesIn.size());
        // Finally, connect to the genome target.
        log.info("Preparing {} genome target {}.", this.targetType, this.outDir);
        this.genomesOut = this.targetType.create(this.outDir, this.clearFlag);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // We count the genomes processed, but full stats for the scores on the changed
        // annotations.
        int gTotal = this.annoMap.size();
        int gCount = 0;
        SummaryStatistics changes = new SummaryStatistics();
        // Loop through the annotation files.
        for (var annoEntry : this.annoMap.entrySet()) {
            gCount++;
            String genomeId = annoEntry.getKey();
            Genome genome = this.genomesIn.getGenome(genomeId);
            if (genome == null)
                throw new IOException("Genome " + genomeId + " not found in " + this.inDir + ".");
            log.info("Processing genome {} of {}: {}.", gCount, gTotal, genome);
            // Set up counters and stats for the genome.
            SummaryStatistics localChanges = new SummaryStatistics();
            int fidCount = 0;
            int skipCount = 0;
            // Loop through the annotation file.
            final File annoFile = annoEntry.getValue();
            try (TabbedLineReader annoStream = new TabbedLineReader(annoFile)) {
                int fidColIdx = annoStream.findField("fid");
                int scoreColIdx = annoStream.findField("score");
                int annoColIdx = annoStream.findField("new_annotation");
                log.info("Reading annotations from file {}.", annoFile);
                for (var line : annoStream) {
                    String fid = line.get(fidColIdx);
                    double score = line.getDouble(scoreColIdx);
                    String annotation = line.get(annoColIdx);
                    fidCount++;
                    // Look for the feature in the genome.
                    Feature feat = genome.getFeature(fid);
                    if (feat == null) {
                        log.error("{} not found in {}.", fid, genome);
                        skipCount++;
                    } else {
                        String old = feat.getPegFunction();
                        if (! annotation.contentEquals(old)) {
                            // Here we have a new annotation.
                            feat.setFunction(annotation);
                            localChanges.addValue(score);
                            changes.addValue(score);
                        }
                    }
                }
                log.info("{} lines read, {} skipped. {} new annotations with mean score {} and score deviation {}.",
                        fidCount, skipCount, localChanges.getN(), localChanges.getMean(),
                        localChanges.getStandardDeviation());
            }
            log.info("Updating genome {}.", genome);
            this.genomesOut.add(genome);
        }
        log.info("{} genomes processed. {} new annotations with mean score {} and score deviation {}.",
                gCount, changes.getN(), changes.getMean(), changes.getStandardDeviation());
    }

}
