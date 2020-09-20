/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.genome.compare.CompareFunctions;
import org.theseed.proteins.Function;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;

/**
 * This command compares genomes in one directory to identically-sequenced genomes in a second directory.  The goal is
 * to create a report of how functional annotations have changed.  An ORF-by-ORF comparison is produced and the
 * functional annotations compared.  A summary of the function mappings is output.
 *
 * The positional parameters are the name of the reference (old) genome directory and the name of the new-genome
 * directory.  The output will display what the old-genome functions are mapped to in the new genomes.
 *
 * The command-line options are as follows.
 *
 * --roles	if specified, a role definition file containing the IDs of important roles
 *
 * @author Bruce Parrello
 *
 */
public class FunctionCompareProcessor extends BaseCompareProcessor {

    /** comparison object */
    private CompareFunctions compareEngine;
    /** optional important-role map */
    private RoleMap roleMap;

    // COMMAND-LINE OPTIONS

    /** important-role definition file */
    @Option(name = "--roles", metaVar = "roles.needed", usage = "important-role definition file")
    private File rolesNeededFile;

    /** new-genome directory */
    @Argument(index = 1, metaVar = "newDir", usage = "new-genome directory")
    private File newDir;

    @Override
    protected void setSubDefaults() {
    }

    /**
     * Validate and load the role map.
     *
     * @throws NoSuchAlgorithmException
     */
    protected void validateSubParms() throws IOException, NoSuchAlgorithmException {
        // Load the role map if needed.
        if (this.rolesNeededFile != null)
            this.roleMap = RoleMap.load(this.rolesNeededFile);
        // Create the comparison engine.
        this.compareEngine = new CompareFunctions();
    }

    @Override
    protected void runCommand() throws Exception {
        // Loop through the new genomes.
        log.info("Scanning new-genome directory {}.", this.newDir);
        GenomeDirectory genomes = new GenomeDirectory(this.newDir);
        log.info("{} genomes found in {}.", genomes.size(), this.newDir);
        for (Genome genome : genomes) {
            File oldGenomeFile = this.findOldGenome(genome);
            if (oldGenomeFile == null)
                log.info("Skipping {}.", genome);
            else {
                log.info("Reading old genome from {} for {}.", oldGenomeFile, genome);
                Genome oldGenome = new Genome(oldGenomeFile);
                boolean found = this.compareEngine.compare(genome, oldGenome);
                if (! found)
                    log.warn("Contig IDs are invalid, comparison for {} and {} aborted.", genome, oldGenome);
            }
        }
        // Now produce the output.
        log.info("Preparing report.");
        System.out.print("old_function\tnew_function\tcount\tpercent");
        if (this.roleMap != null)
            System.out.print("\tneeded");
        System.out.println();
        int pairCount = 0;
        int funCount = 0;
        int simpleCount = 0;
        int emptyCount = 0;
        for (Function oldFun : this.compareEngine.getMissFunctions()) {
            log.debug("Processing {}.", oldFun);
            String funId = oldFun.getId();
            String oldName = oldFun.getName();
            double totalCount = this.compareEngine.getTotalCount(funId);
            // Start this function group.
            int matches = this.compareEngine.getMatchCount(funId);
            System.out.format("%s\t%s\t%d\t%8.2f%n", oldName, "", matches, matches * 100 / totalCount);
            // Loop through the mismatches.
            CountMap<String> missCounts = this.compareEngine.getMissCounts(funId);
            boolean simple = (missCounts.size() == 1 && matches == 0);
            for (CountMap<String>.Count missCount : missCounts.sortedCounts()) {
                String newFun = missCount.getKey();
                matches = missCount.getCount();
                String newName = this.compareEngine.getName(newFun);
                if (newName.isEmpty()) {
                    newName = "(empty string)";
                    simple = false;
                    emptyCount += matches;
                }
                System.out.format("%s\t%s\t%d\t%8.2f", oldName, newName, matches, matches * 100 / totalCount);
                if (this.roleMap != null) {
                    List<Role> roles = Feature.usefulRoles(this.roleMap, newName);
                    String flag = (roles.size() > 0 ? "\tY" : "\t");
                    System.out.print(flag);
                }
                System.out.println();
                pairCount++;
            }
            funCount++;
            if (simple) simpleCount++;
        }
        log.info("{} name pairs found for {} functions. {} functions have simple mappings, {} new-genome proteins had empty functional assignments.",
                pairCount, funCount, simpleCount, emptyCount);
    }

    /**
     * @return the compareEngine
     */
    @Override
    public CompareFunctions getCompareEngine() {
        return this.compareEngine;
    }

}
