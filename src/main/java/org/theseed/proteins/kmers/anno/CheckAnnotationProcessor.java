/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
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
    /** map of confirmed old-to-new annotations */
    private Map<String, String> confirmed;
    /** map of annotation files to process */
    private Map<String, File> annoMap;

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
        // TODO code for validateReporterParms

    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // TODO code for runReporter

    }
}
