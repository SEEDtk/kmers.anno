/**
 *
 */
package org.theseed.proteins.kmers;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.codec.binary.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;

/**
 * This class is used to simplify reading annotation files.  The class itself represents an annotation record,
 * and an iterator is provided to run through the files.
 *
 * Two annotations are the same if the old and new strings match.  The score and feature ID do not count.
 */
public class Annotation {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(Annotation.class);
    /** feature ID */
    private String fid;
    /** score */
    private double score;
    /** old annotation */
    private String oldAnnotation;
    /** new annotation */
    private String newAnnotation;
    /** annotation file name pattern */
    private static final Pattern ANNO_FILE_NAME = Pattern.compile("(\\d+\\.\\d+)\\.anno\\.tbl");

    /**
     * Construct an annotation record.
     *
     * @param fid		feature ID
     * @param scr		annotation score
     * @param oldAnno	old annotation string
     * @param newAnno	new annotation string
     */
    public Annotation(String featureId, double scr, String oldAnno, String newAnno) {
        this.fid = featureId;
        this.score = scr;
        this.newAnnotation = newAnno;
        this.oldAnnotation = oldAnno;
    }

    /**
     * @return TRUE if the annotation is unchanged, else FALSE
     */
    public boolean isGood() {
        return StringUtils.equals(this.newAnnotation, this.oldAnnotation);
    }

    /**
     * @return TRUE if the annotation has a zero score, else FALSE
     */
    public boolean isNull() {
        return Double.isNaN(this.score) || this.score == 0.0;
    }

    /**
     * @return the feature ID
     */
    public String getFid() {
        return this.fid;
    }

    /**
     * @return the score
     */
    public double getScore() {
        return this.score;
    }

    /**
     * @return the old annotation string
     */
    public String getOldAnnotation() {
        return this.oldAnnotation;
    }

    /**
     * @return the new annotation string
     */
    public String getNewAnnotation() {
        return this.newAnnotation;
    }

    /**
     * This is an iterator through an annotation file.
     */
    public static class Iter implements Iterator<Annotation> {

        /** input stream */
        private TabbedLineReader inStream;
        /** feature ID column index */
        private int fidIdx;
        /** score column index */
        private int scoreIdx;
        /** old-annotation column index */
        private int oldAnnoIdx;
        /** new-annotation column index */
        private int newAnnoIdx;

        /**
         * Construct an iterator from an input stream.
         *
         * @param stream	input tabbed line reader
         *
         * @throws IOException
         */
        public Iter(TabbedLineReader stream) throws IOException {
            this.inStream = stream;
            this.fidIdx = stream.findField("fid");
            this.scoreIdx = stream.findField("score");
            this.newAnnoIdx = stream.findField("new_annotation");
            this.oldAnnoIdx = stream.findField("old_annotation");
        }

        @Override
        public boolean hasNext() {
            return this.inStream.hasNext();
        }

        @Override
        public Annotation next() {
            var line = this.inStream.next();
            String fid = line.get(this.fidIdx);
            double score = line.getDouble(this.scoreIdx);
            String newAnno = line.get(this.newAnnoIdx);
            String oldAnno = line.get(this.oldAnnoIdx);
            return new Annotation(fid, score, oldAnno, newAnno);
        }

    }

    /**
     * This is a static method that returns a map of genome IDs to files for an annotation
     * directory.
     *
     * @param annoDir	directory containing the annotation file
     *
     * @return a map from genome IDs to annotation file names
     *
     * @throws IOException
     */
    public static Map<String, File> getAnnoMap(File annoDir) throws IOException {
        // Verify the annotations directory exists.
        if (! annoDir.isDirectory())
            throw new FileNotFoundException("Annotation directory " + annoDir + " is not found or invalid.");
        File[] annoFiles = annoDir.listFiles();
        // Build a map of the annotation files.
        Map<String, File> retVal = new TreeMap<String, File>();
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
                    retVal.put(m.group(1), file);
                }
            }
        }
        log.info("{} annotation files found in {}.", retVal.size(), annoDir);
        return retVal;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.newAnnotation == null) ? 0 : this.newAnnotation.hashCode());
        result = prime * result + ((this.oldAnnotation == null) ? 0 : this.oldAnnotation.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        Annotation other = (Annotation) obj;
        if (this.newAnnotation == null) {
            if (other.newAnnotation != null)
                return false;
        } else if (!this.newAnnotation.equals(other.newAnnotation))
            return false;
        if (this.oldAnnotation == null) {
            if (other.oldAnnotation != null)
                return false;
        } else if (!this.oldAnnotation.equals(other.oldAnnotation))
            return false;
        return true;
    }
}

