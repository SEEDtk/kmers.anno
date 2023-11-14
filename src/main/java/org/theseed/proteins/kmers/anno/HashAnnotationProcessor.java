/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;

import org.theseed.utils.BaseInputProcessor;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.FieldInputStream;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Connection;
import org.theseed.proteins.kmer.hash.ProteinKmerHashMap;

/**
 * This sub-command creates a protein kmer hash and uses it to annotate genomes in a PATRIC
 * genome dump directory.  Such a directory contains a "genome_feature.json" file that has
 * obsolete annotations in it.  The corrected annotations are written to a tab-delimited
 * file of feature IDs and annotations named "anno.tbl".  If no annotation is possible,
 * the old annotation is retained.
 *
 * The standard input should contain a protein annotation file.  This is a tab-delimited file
 * with headers, and the proteins are in a column named "protein", the annotations in a column
 * named "annotation".
 *
 * The positional parameter should be the name of the PATRIC genome dump directory.  This is a
 * directory of sub-directories, each directory containing JSON dump files for one genome.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file containing proteins and annotations (if not STDIN)
 * -K	protein kmer size
 * -b	batch size for load
 *
 * --min	minimum similarity score for an annotation to be acceptable
 * --num	estimated number of proteins being input
 *
 */
public class HashAnnotationProcessor extends BaseInputProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HashAnnotationProcessor.class);
    /** protein map for annotations */
    private ProteinKmerHashMap<String> annoMap;
    /** list of genome directories */
    private File[] gDirs;
    /** input column index for protein strings */
    private int protIdx;
    /** input column index for annotation strings */
    private int annoIdx;

    // COMMAND-LINE OPTIONS

    /** protein kmer size */
    @Option(name = "--kmer", aliases = { "-K" }, metaVar = "10", usage = "protein kmer size")
    private int kmerSize;

    /** batch size for load */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "100", usage = "batch size for hash loading")
    private int batchSize;

    /** minimum similarity score for an annotation */
    @Option(name = "--min", metaVar = "0.1", usage = "minimum acceptable similarity score for annotation")
    private double minScore;

    /** estimated number of proteins being hashed */
    @Option(name = "--num", metaVar = "1000", usage = "estimated number of proteins being put into the hash")
    private int protExpected;

    /** input PATRIC genome dump directory */
    @Argument(index = 0, metaVar = "inDir", usage = "input genome dump directory", required = true)
    private File inDir;

    @Override
    protected void setReaderDefaults() {
        this.kmerSize = 8;
        this.minScore = 0.005;
        this.batchSize = 5000;
        this.protExpected = 15000;
    }

    @Override
    protected void validateReaderParms() throws IOException, ParseFailureException {
        // Insure the kmer size is valid.
        if (this.kmerSize < 2)
            throw new ParseFailureException("Kmer Size must be at least 2.");
        log.info("Kmer size is {}.", this.kmerSize);
        // Insure the batch size is valid.
        if (this.batchSize < 1)
            throw new ParseFailureException("Batch size must be positive.");
        // Insure the estimated protein count is valid.
        if (this.protExpected < 10)
            throw new ParseFailureException("Expected protein count must be at least 10.");
        // Insure the min score is valid.
        if (this.minScore < 0.0 || this.minScore >= 1.0)
            throw new ParseFailureException("Minimum similarity score must be between 0 and 1.");
        log.info("Annotation threshold is {}.", this.minScore);
        // Insure the genome dump directory is valid.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        this.gDirs = this.inDir.listFiles(P3Connection.GENOME_FILTER);
        if (this.gDirs.length <= 0)
            throw new FileNotFoundException("No genomes found in input directory " + this.inDir + ".");
    }

    @Override
    protected void validateReaderInput(TabbedLineReader reader) throws IOException {
        this.protIdx = reader.findField("protein");
        this.annoIdx = reader.findField("annotation");
    }

    @Override
    protected void runReader(TabbedLineReader reader) throws Exception {
        log.info("Creating kmer hash table.");
        this.annoMap = new ProteinKmerHashMap<String>(this.kmerSize, this.protExpected);
        // Initialize the load counters.
        int protCount = 0;
        int skipCount = 0;
        int lineCount = 0;
        // Create the load batch.
        var batch = new ArrayList<Map.Entry<String, String>>(this.batchSize);
        // Loop through the input, storing the proteins.
        for (var line : reader) {
            lineCount++;
            String prot = line.get(this.protIdx);
            String anno = line.get(this.annoIdx);
            if (StringUtils.isBlank(anno) || StringUtils.isBlank(prot)
                    || prot.length() < this.kmerSize || anno.contains("hypothetical"))
                skipCount++;
            else {
                if (batch.size() >= this.batchSize) {
                    // Here the batch is full, so load it.
                    this.loadBatch(batch);
                    batch.clear();
                    log.info("{} lines read, {} proteins kept, {} kmers in table.", lineCount, protCount,
                            this.annoMap.getKmerCount());
                }
                batch.add(new AbstractMap.SimpleEntry<String, String>(prot, anno));
                protCount++;
            }
        }
        // Process the residual batch.
        if (batch.size() > 0)
            this.loadBatch(batch);
        log.info("{} lines read.  {} proteins kept, {} skipped. {} kmers in table.", lineCount,
                protCount, skipCount, this.annoMap.getKmerCount());
        // Now we loop through the genome directories.
        int gCount = 0;
        int foundCount = 0;
        int featCount = 0;
        protCount = 0;
        skipCount = 0;
        long lastMsg = System.currentTimeMillis();
        log.info("Processing {} genome directories.", this.gDirs.length);
        for (File gDir : this.gDirs) {
            String genomeId = gDir.getName();
            gCount++;
            log.info("Processing genome {}:  {}.", gCount, genomeId);
            File featFile = new File(gDir, P3Connection.JSON_FILE_NAME);
            File annoFile = new File(gDir, "anno.tbl");
            try (
                    PrintWriter writer = new PrintWriter(annoFile);
                    FieldInputStream featStream = FieldInputStream.create(featFile)) {
                int fidColIdx = featStream.findField("patric_id");
                int protColIdx = featStream.findField("aa_sequence_md5");
                int annoColIdx = featStream.findField("product");
                // Write the output header.
                writer.println("fid\tannotation");
                // Loop through the genome features.
                for (var line : featStream) {
                    String fid = line.get(fidColIdx);
                    if (StringUtils.isBlank(fid))
                        skipCount++;
                    else {
                        String annotation = line.get(annoColIdx);
                        featCount++;
                        String prot = line.get(protColIdx);
                        if (! StringUtils.isBlank(prot)) {
                            // Here we have a protein to check.
                            protCount++;
                            var closest = this.annoMap.findClosest(prot);
                            if (closest.getSimValue() >= this.minScore) {
                                foundCount++;
                                annotation = closest.getValue();
                            }
                        }
                        writer.println(fid + "\t" + annotation);
                    }
                    if (log.isInfoEnabled() && System.currentTimeMillis() - lastMsg >= 10000) {
                        log.info("{} features processed, [} skipped, {} proteins, {} new annotations found.",
                                featCount, skipCount, protCount, foundCount);
                        lastMsg = System.currentTimeMillis();
                    }
                }
            }
        }
        log.info("{} features processed, [} skipped, {} proteins, {} new annotations found.",
                featCount, skipCount, protCount, foundCount);
    }

    private void loadBatch(ArrayList<Entry<String, String>> batch) {
        batch.parallelStream().forEach(x -> this.annoMap.addProtein(x.getKey(), x.getValue()));
    }

}
