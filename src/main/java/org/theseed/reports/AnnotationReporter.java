/**
 *
 */
package org.theseed.reports;

import java.io.PrintWriter;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;

/**
 * This is the base class for all annotation-comparison reports.  Each such report
 * gets called to initialize and finish, and also once with every old and new version
 * of a feature.
 */
public abstract class AnnotationReporter {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(AnnotationReporter.class);
    /** output print writer */
    private PrintWriter writer;
    /** number of lines written */
    private int counter;
    /** number of fields per line */
    private int width;
    /** output buffer */
    private StringBuffer buffer;

    /**
     * This interface must be supported by the controlling command processor.  The subclass
     * will use it to extract additional parameters.
     */
    public interface IParms {

    }

    /**
     * This enum specifies the type of report.
     */
    public static enum Type {

        FULL {
            @Override
            public AnnotationReporter create() {
                return new FullCompareAnnotationReporter();
            }
        }, NEW_ROLES {
            @Override
            public AnnotationReporter create() {
                return new NewRoleAnnotationReporter();
            }
        };

        /**
         * Create the reporting object for a report of this type.
         */
        public abstract AnnotationReporter create();
    }

    /**
     * Construct an annotation reporter.
     *
     * @param writer	output print writer
     */
    public AnnotationReporter() {
        this.writer = null;
        this.counter = 0;
        this.buffer = new StringBuffer(100);
        this.width = 0;
    }

    /**
     * Write the header line, consisting of multiple fields, tab-delimited.
     *
     * @param fields	array of headers to write
     */
    public void writeHeader(String... fields) {
        String line = StringUtils.join(fields, '\t');
        this.writer.println(line);
        this.width = fields.length;
    }

    /**
     * Write a data line, consisting of multiple fields, tab-delimited.
     *
     * @param fields	array of headers to write
     */
    public synchronized void writeData(Object... fields) {
        this.buffer.setLength(0);
        for (int i = 0; i < this.width; i++) {
            if (buffer.length() >  0)
                buffer.append('\t');
            if (i > fields.length || fields[i] == null)
                this.buffer.append("");
            else
                this.buffer.append(fields[i]);
        }
        this.writer.println(this.buffer.toString());
        this.counter++;
    }

    /**
     * Start the report.
     *
     * @param processor		controlling command processor
     * @param writer		output print writer for the report
     */
    public void startReport(IParms processor, PrintWriter writer) {
        this.writer = writer;
        this.start(processor);
        if (this.width == 0)
            throw new IllegalStateException("AnnotationReporter subclass did not call \"WriteHeader\".");
    }

    /**
     * Start the report.  The subclass must write the header here.
     *
     * @param processor		controlling command processor
     */

    protected abstract void start(IParms processor);

    /**
     * Process an annotation comparison.  Each call here is a data point.
     *
     * @param oldFeat	old version of feature
     * @param newFeat	new version of same feature
     */
    public abstract void processFeature(Feature oldFeat, Feature newFeat);

    /**
     * Finish the report.
     */
    public void finishReport() {
        // Let the subclass do cleanup.
        this.finish();
        // Log the line count.
        log.info("{} lines written to report.", this.counter);
    }

    /**
     * Finish the report.  The subclass should do any needed cleanup here.
     */
    protected abstract void finish();

}
