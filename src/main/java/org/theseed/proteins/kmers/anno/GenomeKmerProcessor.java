/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import org.kohsuke.args4j.Option;
import org.theseed.basic.ICommand;
import org.theseed.genome.Genome;

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
 * -e	minimum number of kmers required for a proposal to be acceptable
 * -i	name of the input file (overrides STDIN)
 * -o	name of the output file (overrides STDOUT)
 *
 * --cache		name of a directory where GTOs are cached
 * --minLength	maximum factor for length decrease when forming a proposal
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

    @Override
    public void setCommandDefaults() {
        this.inFile = null;
        this.outFile = null;
        this.help = false;
    }

    public void runCommand() throws Exception {
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
            newGenome.save(outFile);
        } else {
            log.info("Writing genome to standard output.");
            newGenome.save(System.out);
        }
    }

    @Override
    protected void validateCommandParms() {
    }

}
