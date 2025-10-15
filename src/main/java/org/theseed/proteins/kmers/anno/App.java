package org.theseed.proteins.kmers.anno;


import java.util.Arrays;

import org.theseed.basic.BaseProcessor;
import org.theseed.genome.compare.GeneCopyProcessor;

/**
 * Valid commands are
 *
 *	kmers		annotate a genome using kmer comparison
 *	batch		annotate multiple genomes using kmer comparison
 *	build		build a discriminating-kmer database for a specified list of roles
 *	apply		apply a discriminating-kmer database to genomes to create a role-count file
 *	merge		merge the testing set and the training set into a single file
 *	funMap		map functions between genomes annotated using an old system and newly-annotated genomes
 *	funApply	apply a function mapping to one or more genomes
 *	compare		compare functional assignments between new and old genomes
 *	seqCheck	verify that proteins in genomes are consistently annotated
 *	genes		copy gene names from one genome to a close genome without gene names
 *	hashAnno	use a protein kmer hash to annotate features in a PATRIC dump directory
 *	applyAnno	apply annotations produced by the hash annotator
 *	checkAnno	examine hash-annotator results and write statistics
 *	listAnno	list annotation changes between identical genomes
 *	updateJson	update annotations in JSON genome files
 *	buildGtos	build GTOs from PATRIC data and annotation update files
 */
public class App
{

    protected static final String[] COMMANDS = new String[] {
             "kmers", "annotate a genome using kmer comparison",
             "batch", "annotate multiple genomes using kmer comparison",
             "build", "build a discriminating-kmer database for a specified list of roles",
             "apply", "apply a discriminating-kmer database to genomes to create a role-count file",
             "merge", "merge the testing set and the training set into a single file",
             "funMap", "map functions between genomes annotated using an old system and newly-annotated genomes",
             "funApply", "apply a function mapping to one or more genomes",
             "compare", "compare functional assignments between new and old genomes",
             "seqCheck", "verify that proteins in genomes are consistently annotated",
             "genes", "copy gene names from one genome to a close genome without gene names",
             "hashAnno", "use a protein kmer hash to annotate features in a PATRIC dump directory",
             "applyAnno", "apply annotations produced by the hash annotator",
             "checkAnno", "examine hash-annotator results and write statistics",
             "listAnno", "list annotation changes between identical genomes",
             "updateJson", "update annotations in JSON genome files",
             "buildGtos", "build GTOs from PATRIC data and annotation update files"
    };

    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        // Parse the parameters.
        switch (command) {
        case "kmers" -> processor = new GenomeKmerProcessor();
        case "batch" -> processor = new BatchKmerProcessor();
        case "build" -> processor = new BuildKmerProcessor();
        case "apply" -> processor = new ApplyKmerProcessor();
        case "merge" -> processor = new MergeFilesProcessor();
        case "funMap" -> processor = new FunctionCompareProcessor();
        case "funApply" -> processor = new FunctionApplyProcessor();
        case "compare" -> processor = new GenomeCompareProcessor();
        case "seqCheck" -> processor = new SequenceCheckProcessor();
        case "genes" -> processor = new GeneCopyProcessor();
        case "hashAnno" -> processor = new HashAnnotationProcessor();
        case "applyAnno" -> processor = new ApplyAnnotationProcessor();
        case "checkAnno" -> processor = new CheckAnnotationProcessor();
        case "listAnno" -> processor = new ListNewAnnotationProcessor();
        case "updateJson" -> processor = new UpdateJsonProcessor();
        case "buildGtos" -> processor = new GtoBuildProcessor();
        case "-h", "--help" -> processor = null;
        default -> throw new RuntimeException("Invalid command " + command + ".");
        }
        if (processor == null)
            BaseProcessor.showCommands(COMMANDS);
        else {
            processor.parseCommand(newArgs);
            processor.run();
        }
    }
}
