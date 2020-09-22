/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.genome.compare.CompareType;
import org.theseed.genome.compare.IGenomeMatcher;
import org.theseed.genome.compare.MatchGenomes;

/**
 * This command performs a functional assignment comparison between identically-sequenced genomes.  A directory of reference genomes is specified
 * along with one or more directories of newly-annotated genomes.  For each newly-annotated genome, we find a sequence-identical reference genome
 * and determine the percentage of functional assignment matches.  The output is a matrix, with each column representing a new-genome directory
 * and each row a corresponding reference genome.
 *
 * The positional parameters are the name of the reference-genome directory followed by the names of the new-genome directories.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more detailed progress messages
 * -t 	type of comparison (default FUNCTIONS)
 *
 * @author Bruce Parrello
 *
 */
public class GenomeCompareProcessor extends BaseCompareProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeCompareProcessor.class);
    /** comparison engine */
    private IGenomeMatcher compareEngine;
    /** map of old-genome IDs to quality percentages, corresponding to the input directory positions */
    private Map<String, String[]> genomeMatchMap;
    /** number of good matches in each directory */
    private int good[];
    /** number of bad matches in each directory */
    private int bad[];

    // COMMAND-LINE OPTIONS

    /** type of comparison */
    @Option(name = "-t", aliases = { "--type" }, usage = "type of comparison to perform")
    private CompareType type;

    /** new-genome directories */
    @Argument(index = 1, metaVar = "newDir1 newDir2 ...", usage = "directory of new (modified) genomes", multiValued = true)
    private List<File> newDirs;

    @Override
    protected void setSubDefaults() {
    }

    @Override
    protected MatchGenomes getCompareEngine() {
        return (MatchGenomes) this.compareEngine;
    }

    @Override
    protected void validateSubParms() throws IOException, NoSuchAlgorithmException {
        // Create the comparison engine.
        this.compareEngine = this.type.create();
        // Verify the new-genome directories.
        for (File newDir : this.newDirs) {
            if (! newDir.isDirectory())
                throw new FileNotFoundException("New-genome directory " + newDir + " is not found or invalid.");
        }
        // Create the genome-matching map.
        this.genomeMatchMap = new TreeMap<String, String[]>();
        // Initialize the totals.
        this.good = new int[this.newDirs.size()];
        this.bad = new int[this.newDirs.size()];
    }

    @Override
    protected void runCommand() throws Exception {
        // This will track our position in the directories array.
        int iDir = 0;
        // We will loop through the new-genome directories, processing each one individually.
        for (File newDir : this.newDirs) {
            log.info("Processing input directory {}.", newDir);
            GenomeDirectory genomes = new GenomeDirectory(newDir);
            for (Genome genome : genomes) {
                // Locate the old genome for this new one.
                File oldFile = this.findOldGenome(genome);
                if (oldFile == null)
                    log.warn("No reference match for {}-- skipping.", genome);
                else {
                    Genome oldGenome = new Genome(oldFile);
                    // Note the old genome goes first in the check call, since we may end up updating the second
                    // one, and we don't want to update the old genome, which is reused.
                    log.info("Comparing {} to {}.", genome, oldGenome);
                    boolean ok = this.compareEngine.compare(oldGenome, genome);
                    if (! ok)
                        log.error("Contig IDs in {} are invalid.  Comparison aborted.", genome);
                    else {
                        // Here we have a result.
                        String[] resultArray = this.genomeMatchMap.computeIfAbsent(oldGenome.getId(), x -> new String[this.newDirs.size()]);
                        resultArray[iDir] = String.format("%8.4f", this.compareEngine.percent());
                        this.good[iDir] += this.compareEngine.getGood();
                        this.bad[iDir] += this.compareEngine.getBad();
                    }
                }
            }
            iDir++;
        }
        log.info("Writing report.");
        // The tricky part of this report is that nulls have to be converted to empty strings.
        System.out.println(this.newDirs.stream().map(x -> x.getName()).collect(Collectors.joining("\t", "reference\t", "")));
        for (Map.Entry<String, String[]> percentEntry : this.genomeMatchMap.entrySet()) {
            String[] percents = percentEntry.getValue();
            String refGenome = percentEntry.getKey();
            System.out.println(Arrays.stream(percents).map((x -> (x == null ? "" : x))).collect(Collectors.joining("\t", refGenome + "\t", "")));
        }
        // Output the totals.
        System.out.println();
        String[] percents = new String[this.newDirs.size()];
        for (int i = 0; i < percents.length; i++) {
            double pct = 0.0;
            if (this.good[i] > 0) {
                pct = this.good[i] * 100.0 / (this.bad[i] + this.good[i]);
                percents[i] = String.format("%8.4f", pct);
            }
        }
        System.out.println(Arrays.stream(percents).collect(Collectors.joining("\t", "TOTAL\t", "")));
    }

}
