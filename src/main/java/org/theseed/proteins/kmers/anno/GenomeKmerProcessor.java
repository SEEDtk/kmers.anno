/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.IOException;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.genome.Genome;
import org.theseed.utils.ICommand;

/**
 * This program attempts to project proteins from close genomes to a new genome.  The new genome should have its taxonomy
 * information, close-genome list, genetic code, ID, and name filled in, as well as all of its contigs.  We loop through
 * the N closest genomes, using kmers to project proteins onto ORFs.
 *
 * The genome should be in the form of a GTO.  It is read from the standard input, and an annotated GTO is written to the
 * standard output.
 *
 * The following command-line options are supported.
 *
 * -m	minimum strength for an annotation proposal to be acceptable
 * -f	maximum factor for length increase in forming a proposal
 * -K	kmer length to use
 * -n	number of close genomes to use
 * -i	name of the input file (overrides STDIN)
 * -o	name of the output file (overrides STDOUT)
 *
 * --algorithm	kmer algorithm to use-- STRICT (only match kmers unique in the contigs) or AGGRESSIVE (match all kmers in the contigs)
 * --trace		if specified, the name of a functional assignment; information about the function will be written to the log
 *
 * @author Bruce Parrello
 *
 */
public class GenomeKmerProcessor extends KmerProcessor implements ICommand {

    // COMMAND LINE

    /** input file name */
    @Option(name = "-i", aliases = { "--input" }, usage = "input file name (if not STDIN)")
    private File inFile;

    /** output file name */
    @Option(name = "-o", aliases = { "--output" }, usage = "output file name (if not STDOUT)")
    private File outFile;

    /**
     * Parse command-line options to specify the parameters of this object.
     *
     * @param args	an array of the command-line parameters and options
     *
     * @return TRUE if successful, FALSE if the parameters are invalid
     */
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.inFile = null;
        this.outFile = null;
        this.help = false;
        setDefaults();
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                validateParms();
                // Denote we can run.
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            // For parameter errors, we display the command usage.
            parser.printUsage(System.err);
        }
        return retVal;
    }

    public void run() {
        try {
            log.info("Connecting to PATRIC.");
            Genome newGenome;
            if (this.inFile != null) {
                log.info("Reading genome from {}.", this.inFile);
                newGenome = new Genome(this.inFile);
            } else {
                log.info("Reading genome from standard input.");
                newGenome = new Genome(System.in);
            }
            annotateGenome(newGenome);
            // Write the result.
            if (this.outFile != null) {
                log.info("Writing genome to {}.", this.outFile);
                newGenome.update(outFile);
            } else {
                log.info("Writing genome to standard output.");
                newGenome.update(System.out);
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
