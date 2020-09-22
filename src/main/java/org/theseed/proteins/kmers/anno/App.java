package org.theseed.proteins.kmers.anno;


import java.util.Arrays;

import org.theseed.utils.BaseProcessor;

/**
 * Valid commands are
 *
 * 	kmers		annotate a genome using kmer comparison
 * 	batch		annotate multiple genomes using kmer comparison
 * 	build		build a discriminating-kmer database for a specified list of roles
 * 	apply		apply a discriminating-kmer database to genomes to create a role-count file
 * 	merge		merge the testing set and the training set into a single file
 * 	funMap		map functions between genomes annotated using an old system and newly-annotated genomes
 *	funApply	apply a function mapping to one or more genomes
 *	compare		compare functional assignments between new and old genomes
 *	subCompare	compare subsystems between two sets of genomes
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        // Parse the parameters.
        switch (command) {
        case "kmers" :
            processor = new GenomeKmerProcessor();
            break;
        case "batch" :
            processor = new BatchKmerProcessor();
            break;
        case "build" :
            processor = new BuildKmerProcessor();
            break;
        case "apply" :
            processor = new ApplyKmerProcessor();
            break;
        case "merge" :
            processor = new MergeFilesProcessor();
            break;
        case "funMap" :
            processor = new FunctionCompareProcessor();
            break;
        case "funApply" :
            processor = new FunctionApplyProcessor();
            break;
        case "compare" :
            processor = new GenomeCompareProcessor();
            break;
        default:
            throw new RuntimeException("Invalid command \"" + command + "\".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
