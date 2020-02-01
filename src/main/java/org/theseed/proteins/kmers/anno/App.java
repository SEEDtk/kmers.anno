package org.theseed.proteins.kmers.anno;

import org.theseed.utils.ICommand;

/**
 * Main application
 *
 */
public class App
{
    public static void main( String[] args )
    {
        ICommand runObject = new KmerProcessor();
        boolean ok = runObject.parseCommand(args);
        if (ok) {
            runObject.run();
        }
    }
}
