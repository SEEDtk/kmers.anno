/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.kmers.KmerReference;
import org.theseed.reports.ApplyKmerReporter;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.BaseProcessor;

/**
 * This command applies a discriminating-kmer database to one or more genomes to count the occurrences of interesting roles.
 * The output file can be used as a testing set for the TrainProcessor in dl4j.eval.
 *
 * The positional parameters are the name of the file containing the discriminating kmers, the name of the roles-in-use file
 * (which determines ordering for the output columns), and the name of the input GTO directory.
 *
 * The command-line parameters are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent progress messages on the log
 * -m	mimimum number of required hits
 *
 * --format		output format (VERIFY or TRAIN)
 *
 * @author Bruce Parrello
 *
 */
public class ApplyKmerProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ApplyKmerProcessor.class);
    /** reporting facility */
    private ApplyKmerReporter reporter;
    /** main kmer hash database */
    private Map<String, String> kmerRoleMap;

    // COMAND-LINE OPTIONS

    /** report output format */
    @Option(name = "--format", usage = "reporting format")
    private ApplyKmerReporter.Type outputType;

    @Option(name = "-m", aliases = { "--min" }, metaVar = "10", usage = "minimum number of hits required to call a role")
    private int minHits;

    /** kmer input file */
    @Argument(index = 0, metaVar = "kmerdb.tbl", usage = " discriminating kmer database")
    private File kmerDbFile;

    /** roles-in-use file */
    @Argument(index = 1, metaVar = "roles.in.use", usage = "list of roles in use")
    private File goodRoleFile;

    /** input GTO directory */
    @Argument(index = 2, metaVar = "gtoDir", usage = "input genome directory")
    private File inDir;

    @Override
    protected void setDefaults() {
        this.outputType = ApplyKmerReporter.Type.APPLY;
        this.minHits = 5;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Verify the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " not found or invalid.");
        // Verify the kmer database.
        if (! this.kmerDbFile.canRead())
            throw new FileNotFoundException("Kmer database file " + this.kmerDbFile + " not found or unreadable.");
        // Verify the minimum number of hits.
        if (this.minHits < 1)
            throw new IllegalArgumentException("Min-hits must be positive.");
        // Initialize the reporting.
        this.reporter = this.outputType.create(System.out);
        if (! this.goodRoleFile.canRead())
            throw new FileNotFoundException("Roles-to-use file " + this.goodRoleFile + " not found or unreadable.");
        log.info("Reading roles to use from {}.", this.goodRoleFile);
        this.reporter.initReport(this.goodRoleFile);
        // Load the kmer database.
        log.info("Loading kmer database from {}.", this.kmerDbFile);
        this.kmerRoleMap = new HashMap<String, String>((int) (this.kmerDbFile.length() / 30));
        try (TabbedLineReader kmerStream = new TabbedLineReader(this.kmerDbFile, 2)) {
            String kmer = "";
            for (TabbedLineReader.Line line : kmerStream) {
                kmer = line.get(0);
                this.kmerRoleMap.put(kmer, line.get(1));
            }
            KmerReference.setKmerSize(kmer.length());
            log.info("Kmer size is {}.", KmerReference.getKmerSize());
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Loop through the genomes.
        GenomeDirectory genomes = new GenomeDirectory(this.inDir);
        log.info("{} genomes found in input directory.", genomes.size());
        for (Genome genome : genomes) {
            log.info("Processing genome {}.", genome);
            this.reporter.openGenome(genome);
            // Loop through the proteins.
            for (Feature feat : genome.getPegs()) {
                ProteinKmers kmers = new ProteinKmers(feat.getProteinTranslation());
                // We consider ourselves to have a hit if the kmers all hit one role.
                String roleId = null;
                int count = 0;
                boolean badPeg = false;
                Iterator<String> iter = kmers.iterator();
                while (iter.hasNext() && ! badPeg) {
                    String possible = this.kmerRoleMap.get(iter.next());
                    if (possible != null) {
                        // Here we have a known kmer.
                        if (roleId == null) {
                            // This is the first hit for the feature.
                            roleId = possible;
                            count = 1;
                        } else if (possible.contentEquals(roleId)) {
                            // This is a confirming hit.
                            count++;
                        } else {
                            // This is a bad hit; the peg is ambiguous.
                            badPeg = true;
                        }
                    }
                }
                if (roleId != null && ! badPeg && count >= this.minHits)
                    this.reporter.recordFeature(feat, roleId, count);
            }
            // Close this genome.
            this.reporter.closeGenome();
        }
        // Close the report.
        this.reporter.closeReport();
        this.reporter.close();
    }

}
