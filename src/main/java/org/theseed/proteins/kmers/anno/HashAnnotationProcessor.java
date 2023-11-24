/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Duration;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.kmer.hash.GenomeProteinKmers;
import org.theseed.utils.BaseMultiReportProcessor;

/**
 * This sub-command uses protein kmer hashing to annotate GTOs.  The GTOs should be at least
 * PROTEIN detail level:  DNA is not required.  For each GTO, a file named "XXXXXX.anno.tbl"
 * will be output, where "XXXXXX" is the genome ID. If no annotation is possible, the old annotation
 * is retained.
 *
 * The annotation data is in a role annotation file.  This is a tab-delimited file with headers, and
 * the proteins are in a column named "protein" and the annotations in a column named "annotation".
 *
 * The positional parameters should be the name of the role annotation file and a genome source.
 *
 * Our strategy is to read in all the kmers for a genome, then process the entire role annotation file
 * to find the best match for each incoming peg.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -K	protein kmer size (default 8)
 * -D	output directory for annotation files (default "Annotations" in the current directory)
 * -t	type of genome source (default DIR)
 *
 * --minSim		minimum similarity score for an annotation to be acceptable (default 0.005)
 * --minLen		minimum acceptable protein length (default 50)
 *
 */
public class HashAnnotationProcessor extends BaseMultiReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HashAnnotationProcessor.class);
    /** input column index for protein strings */
    private int protIdx;
    /** input column index for annotation strings */
    private int annoIdx;
    /** genomes to process */
    private GenomeSource genomes;

    // COMMAND-LINE OPTIONS

    /** protein kmer size */
    @Option(name = "--kmer", aliases = { "-K" }, metaVar = "10", usage = "protein kmer size")
    private int kmerSize;

    /** minimum similarity score for an annotation */
    @Option(name = "--minSim", metaVar = "0.1", usage = "minimum acceptable similarity score for annotation")
    private double minScore;

    /** minimum length for a usable annotation protein */
    @Option(name = "--minLen", metaVar = "200", usage = "minimum acceptable length for an annotation protein")
    private double minProtLen;

    /** type of genome source */
    @Option(name = "--source", aliases = { "-t" }, usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** input role annotation file */
    @Argument(index = 0, metaVar = "annoFile", usage = "input role annotation file", required = true)
    private File annoFile;

    /** input genome source */
    @Argument(index = 1, metaVar = "inDir", usage = "input genome source", required = true)
    private File inDir;

    @Override
    protected void setMultiReportDefaults() {
        this.kmerSize = 8;
        this.minScore = 0.005;
        this.sourceType = GenomeSource.Type.DIR;
        this.minProtLen = 50;
    }
    @Override
    protected File setDefaultOutputDir(File curDir) {
        return new File(curDir, "Annotations");
    }

    @Override
    protected void validateMultiReportParms() throws IOException, ParseFailureException {
        // Insure the kmer size is valid.
        if (this.kmerSize < 2)
            throw new ParseFailureException("Kmer Size must be at least 2.");
        log.info("Kmer size is {}.", this.kmerSize);
        // Insure the min score is valid.
        if (this.minScore < 0.0 || this.minScore >= 1.0)
            throw new ParseFailureException("Minimum similarity score must be between 0 and 1.");
        log.info("Annotation threshold is {}.", this.minScore);
        // Insure the genome dump directory is valid.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        // Insure the minimum protein length is valid.
        if (this.minProtLen < this.kmerSize)
            throw new FileNotFoundException("Minimum protein length cannot be less than kmer size.");
        // Validate the role annotation file.  We make sure we can read it and that it has the necessary fields.
        if (! this.annoFile.canRead())
            throw new FileNotFoundException("Role annotation file " + this.annoFile + " is not found or unreadable.");
        try (TabbedLineReader annoStream = new TabbedLineReader(this.annoFile)) {
            this.annoIdx = annoStream.findField("annotation");
            this.protIdx = annoStream.findField("protein");
        }
        // Validate the genome source.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Genome source " + this.inDir + " not found.");
        log.info("Connecting to genome source {} of type {}.", this.inDir, this.sourceType);
        this.genomes = this.sourceType.create(this.inDir);
        log.info("{} genomes loaded from {}.", this.genomes.size(), this.inDir);
    }

    @Override
    protected void runMultiReports() throws Exception {
        // Now we loop through the genome directories.
        int gTotal = this.genomes.size();
        int gCount = 0;
        int featCount = 0;
        int protCount = 0;
        int confirmCount = 0;
        int defaultCount = 0;
        int newAnnoCount = 0;
        for (Genome genome : this.genomes) {
            String genomeId = genome.getId();
            gCount++;
            log.info("Processing genome {} of {}:  {}.", gCount, gTotal, genome);
            try (PrintWriter writer = this.openReport(genomeId + ".anno.tbl")) {
                long genomeStart = System.currentTimeMillis();
                // Loop through the genome features, filling the kmer hash.
                int fCount = 0;
                int sCount = 0;
                int pCount = 0;
                GenomeProteinKmers genomeKmers = new GenomeProteinKmers(this.kmerSize, this.minScore);
                for (Feature feat : genome.getFeatures()) {
                    // Get the data for this feature.
                    String prot = feat.getProteinTranslation();
                    fCount++;
                    if (StringUtils.isBlank(prot))
                        sCount++;
                    else {
                        String fid = feat.getId();
                        String annotation = feat.getPegFunction();
                        pCount++;
                        genomeKmers.addProtein(fid, prot, annotation);
                    }
                }
                featCount += fCount;
                protCount += pCount;
                log.info("{} kmers for {} proteins in hash for {}.", genomeKmers.getKmerCount(), pCount, genome);
                log.info("{} features processed, {} skipped, {} proteins", fCount, sCount, pCount);
                // Now we need to annotate the features.  We do this by reading through the annotation file.
                try (TabbedLineReader annoStream = new TabbedLineReader(this.annoFile)) {
                    int closeCount = 0;
                    int lineCount = 0;
                    long lastMsg = System.currentTimeMillis();
                    long start = lastMsg;
                    // Loop through the annotation stream, performing annotations.
                    for (var line : annoStream) {
                        String prot = line.get(this.protIdx);
                        String anno = line.get(this.annoIdx);
                        lineCount++;
                        // Insure this is a useful annotation.
                        if (! StringUtils.isBlank(anno) && prot.length() >= this.minProtLen) {
                            int count = genomeKmers.processProposal(prot, anno);
                            closeCount += count;
                            if (log.isInfoEnabled() && System.currentTimeMillis() - lastMsg >= 5000) {
                                lastMsg = System.currentTimeMillis();
                                double rate = lineCount * 1000.0 / (lastMsg - start);
                                log.info("{} annotation lines read, {} matches found, {} lines/second.", lineCount, closeCount, rate);
                            }
                        }
                    }
                    // Write the output header.
                    writer.println("fid\tscore\tnew_annotation\told_annotation");
                    // Now unspool the features to the output.  Note that we will leave the score blank if we
                    // have a non-protein feature.
                    int cCount = 0;
                    int naCount = 0;
                    int dCount = 0;
                    for (Feature feat : genome.getFeatures()) {
                        String md5 = feat.getMD5();
                        String oldAnnotation = feat.getPegFunction();
                        GenomeProteinKmers.Proposal proposal = null;
                        if (! StringUtils.isBlank(md5))
                            proposal = genomeKmers.getProposal(md5);
                        // At this point, a NULL proposal means the feature was not annotated.
                        if (proposal == null)
                            writer.println(feat.getId() + "\t\t" + oldAnnotation + "\t" + oldAnnotation);
                        else {
                            // Here we have a successful annotation.  Determine the type (default, confirm, new).
                            double score = proposal.getSim();
                            String newAnnotation = proposal.getAnnotation();
                            if (score == 0.0)
                                dCount++;
                            else if (oldAnnotation.contentEquals(newAnnotation))
                                cCount++;
                            else
                                naCount++;
                            writer.println(feat.getId() + "\t" + String.valueOf(score) + "\t" + newAnnotation
                                    + "\t" + oldAnnotation);
                        }
                    }
                    if (log.isInfoEnabled()) {
                        Duration time = Duration.ofMillis(System.currentTimeMillis() - genomeStart);
                        log.info("{} default annotations, {} confirmed annotations, {} new annotations in {}.",
                                dCount, cCount, naCount, genome);
                        log.info("{} to annotation {}.", time, genome);
                    }
                    confirmCount += cCount;
                    defaultCount += dCount;
                    newAnnoCount += naCount;
                }
            }
            log.info("{} total proteins out of {} features processed for {} genomes.", protCount, featCount, gCount);
            log.info("{} annotations confirmed, {} updated, {} defaulted.", confirmCount, newAnnoCount, defaultCount);
        }
    }

}
