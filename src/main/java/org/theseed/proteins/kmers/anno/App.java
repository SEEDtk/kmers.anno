package org.theseed.proteins.kmers.anno;


/**
 * Hello world!
 *
 */
public class App
{
    public static void main( String[] args )
    {
        Processor runObject = new Processor();
        boolean ok = runObject.parseCommand(args);
        if (ok) {
            runObject.run();
        }
    }
}
