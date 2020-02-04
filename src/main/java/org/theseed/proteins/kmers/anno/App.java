package org.theseed.proteins.kmers.anno;


import java.util.Arrays;

import org.theseed.utils.ICommand;

/**
 * Valid commands are
 *
 * 		kmers		annotate a genome using kmer comparison
 * 		compare		compare the called ORFs of two GTOs
 *
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        ICommand processor;
        // Parse the parameters.
        switch (command) {
        case "annotate" :
            processor = new GenomeKmerProcessor();
            break;
        case "rerun" :
            processor = new BatchKmerProcessor();
            break;
        default:
            throw new RuntimeException("Invalid command.  Must be \"kmers\".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
