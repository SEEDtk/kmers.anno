package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.ICommand;

/**
 * This program attempts to project proteins from close genomes to genomes that need re-annotation.  For each one,
 * we loop through  * the N closest genomes, using kmers to project proteins onto ORFs.
 *
 * The incoming genomes should be in the form of GTOs.  We read a list of input and output file names from
 * an input file.  These come in via a tab-delimited file, with the input names in the first column and the
 * output names in the second.  The file names should be relative to the input file's directory.
 *
 * The positional parameter is the name of the input file.
 *
 * The following command-line options are supported.
 *
 * -m	minimum strength for an annotation proposal to be acceptable
 * -f	maximum factor for length increase in forming a proposal
 * -K	kmer length to use
 * -n	number of close genomes to use
 * -e	minimum number of kmers required for a proposal to be acceptable
 *
 * --cache		name of a directory where GTOs are cached
 * --minLength	maximum factor for length decrease when forming a proposal
 * --algorithm	kmer algorithm to use-- STRICT (only match kmers unique in the contigs) or AGGRESSIVE (match all kmers in the contigs)
 * --trace		if specified, the name of a functional assignment; information about the function will be written to the log
 *
 * @author Bruce Parrello
 *
 */
public class BatchKmerProcessor extends KmerProcessor implements ICommand {

    // COMMAND-LINE OPTIONS

    @Argument(index = 0, usage = "input file containing input and output GTO names", required = true)
    private File inFile;

    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        setDefaults();
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                validateParms();
                // Verify the input file.
                if (! this.inFile.canRead())
                    throw new IOException("Input file " + this.inFile + " not found or unreadable.");
                // Denote we can run.
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            // For parameter errors, we display the command usage.
            parser.printUsage(System.err);
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
        return retVal;
    }

    @Override
    public void run() {
        try {
            // The basic process is to run through the input file one line at a time.
            long start = System.currentTimeMillis();
            File dir = this.inFile.getAbsoluteFile().getParentFile();
            log.info("Reading GTO names from {} in directory {}.", this.inFile, dir);
            int gCount = 0;
            try (TabbedLineReader reader = new TabbedLineReader(this.inFile, 2)) {
                // Loop through the input.
                for (TabbedLineReader.Line line : reader) {
                    // The output file name is in the second column.
                    File outFile = new File(dir, line.get(1));
                    // Read the genome from the first column.
                    File genomeFile = new File(dir, line.get(0));
                    log.info("Reading genome from {}.", genomeFile);
                    Genome genome = new Genome(genomeFile);
                    // Remove all the current annotations.
                    genome.deAnnotate();
                    // Annotate the genome.
                    this.annotateGenome(genome);
                    // Write it to the output.
                    log.info("Writing genome to {}.", outFile);
                    genome.update(outFile);
                    gCount++;
                }
            }
            log.info("Processing complete.  {} genomes annotated, {} seconds / genome.", gCount, (double) (System.currentTimeMillis() - start) / (gCount * 1000));
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
