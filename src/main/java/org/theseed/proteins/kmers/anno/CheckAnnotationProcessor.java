/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.kmers.Annotation;
import org.theseed.utils.BaseReportProcessor;

/**
 * This sub-command analyzes the results of a hash-annotation run.  The hash-annotator produces an XXXXXX.X.anno.tbl file
 * for each annotated genome XXXXXX.X, along with a "changes.tbl" file that contains all the annotations that changed.
 * This program reads the "changes.tbl" file and memorizes confirmed re-annotations, that is, annotation mappings that
 * have very high scores.  For each output genome, we then split the annotation recommendations into those where the
 * new annotation has a score of zero, those where the new annotation is the same as the old one or is confirmed, and
 * those where the new annotation is unconfirmed.  We then produce statistics on the scores for the latter two groups
 * for each genome.  These will be displayed on the standard output.
 *
 * The positional parameter is the name of the annotation directory.  The command-line options are as follows:
 *
 *  -h	display command-line usage
 *  -v	display more frequent log messages
 *  -o	output file for report (if not STDOUT)
 *  -m	minimum score for a confirmed annotation change (default 0.9)
 *
 * @author Bruce Parrello
 *
 */
public class CheckAnnotationProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(CheckAnnotationProcessor.class);
    /** map of confirmed new-to-old annotations */
    private Set<Annotation> confirmed;
    /** map of annotation files to process */
    private Map<String, File> annoMap;
    /** statistics for good changes */
    private SummaryStatistics goodStats;
    /** statistics for questionable changes */
    private SummaryStatistics badStats;
    /** count of null changes */
    private int keepCount;
    private int featCount;

    // COMMAND-LINE OPTIONS

    /** minimum confirmed-change score */
    @Option(name = "--min", aliases = { "-m" }, metaVar = "0.95", usage = "minimum score for a confirmed re-annotation")
    private double minScore;

    /** input annotation directory */
    @Argument(index = 0, metaVar = "annoDir", usage = "input annotation directory")
    private File inDir;

    @Override
    protected void setReporterDefaults() {
        this.minScore = 0.9;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify the minimum score.
        if (this.minScore <= 0.0 || this.minScore > 1.0)
            throw new ParseFailureException("Minimum score must be greater than 0 and no greater than 1.");
        // Set up the input directory.
        this.annoMap = Annotation.getAnnoMap(this.inDir);
        // Insure we can read the changes.
        File changeFile = new File(this.inDir, "changes.tbl");
        if (! changeFile.canRead())
            throw new FileNotFoundException("Annotation directory " + this.inDir + " does not have a readable changes file.");
        // Now create the confirm list.  This can be construed as a set, because the key of an annotation
        // is the two annotation strings.
        this.confirmed = new HashSet<Annotation>(4000);
        log.info("Reading confirmed changes from {}.", changeFile);
        try (TabbedLineReader changeStream = new TabbedLineReader(changeFile)) {
            Iterator<Annotation> iter = new Annotation.Iter(changeStream);
            int changeCount = 0;
            while (iter.hasNext()) {
                Annotation anno = iter.next();
                changeCount++;
                // If the score is high enough for a confirmed annotation, so we save it.
                if (anno.getScore() >= this.minScore)
                    this.confirmed.add(anno);
            }
            log.info("{} changes checked, {} were confirmed.", changeCount, this.confirmed.size());
        }
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Initialize the global counters and stats.
        this.keepCount = 0;
        this.featCount = 0;
        this.goodStats = new SummaryStatistics();
        this.badStats = new SummaryStatistics();
        // Write the output header.
        writer.println("genome\tfids\tdefaulted\tgood_count\tgood_mean\tgood_min\tgood_sdev\tother_count\tother_mean\tother_min\tother_sdev");
        // Loop through the genomes, processing the annotation files.
        int gCount = 0;
        final int gTotal = this.annoMap.size();
        for (var annoEntry : this.annoMap.entrySet()) {
            gCount++;
            String genomeId = annoEntry.getKey();
            File annoFile = annoEntry.getValue();
            log.info("Processing genome {} of {}:  {} from {}.", gCount, gTotal, genomeId, annoFile);
            // Set up the genome-level statistics.
            SummaryStatistics good = new SummaryStatistics();
            SummaryStatistics bad = new SummaryStatistics();
            int keep = 0;
            int feat = 0;
            // Now open and process the file.
            try (TabbedLineReader annoStream = new TabbedLineReader(annoFile)) {
                Iterator<Annotation> iter = new Annotation.Iter(annoStream);
                while (iter.hasNext()) {
                    feat++;
                    this.featCount++;
                    Annotation anno = iter.next();
                    double score = anno.getScore();
                    if (anno.isNull()) {
                        this.keepCount++;
                        keep++;
                    } else if (anno.isGood() || this.confirmed.contains(anno)) {
                        good.addValue(score);
                        this.goodStats.addValue(score);
                    } else {
                        bad.addValue(score);
                        this.badStats.addValue(score);
                    }
                }
                this.report(writer, genomeId, feat, keep, good, bad);
            }
        }
        // Write the totals.
        this.report(writer, "TOTALS", this.featCount, this.keepCount, this.goodStats, this.badStats);
        log.info("{} genomes processed.", gCount);
    }

    /**
     * Write a line of statistics to the output.
     *
     * @param writer	output print writer
     * @param genomeId	ID of the genome
     * @param feat		number of features processed
     * @param keep		number of default annotations
     * @param good		summary statistics for known good annotations
     * @param bad		summary statistics for questionable annotations
     */
    private void report(PrintWriter writer, String genomeId, int feat, int keep, SummaryStatistics good,
            SummaryStatistics bad) {
        writer.println(genomeId + "\t" + Integer.toString(feat) + "\t" + Integer.toString(keep)
                + "\t" + Long.toString(good.getN()) + "\t" + Double.toString(good.getMean())
                + "\t" + Double.toString(good.getMin()) + "\t" + Double.toString(good.getStandardDeviation())
                + "\t" + Long.toString(bad.getN()) + "\t" + Double.toString(bad.getMean())
                + "\t" + Double.toString(bad.getMin()) + "\t" + Double.toString(bad.getStandardDeviation())
                );
    }
}
