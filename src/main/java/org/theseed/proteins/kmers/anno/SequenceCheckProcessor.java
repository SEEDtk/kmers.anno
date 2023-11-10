/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.MD5Hex;

/**
 * This command determines whether or not proteins in a set of genomes are consistently annotated.  Each genome's protein features will
 * be organized by protein sequence MD5.  If two features with the same sequence have different annotations, the difference is flagged.
 *
 * This process requires a staggering amount of memory.
 *
 * The positional parameter is the name of the directory containing the genomes.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more detailed progress messages
 *
 * --roles	role definition file containing roles in subsystems
 *
 * @author Bruce Parrello
 *
 */
public class SequenceCheckProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SequenceCheckProcessor.class);
    /** map from protein ID to feature list */
    private Map<String, Set<Feature>> proteinMap;
    /** interesting role map */
    private RoleMap roleMap;

    // COMMAND-LINE OPTIONS

    @Option(name = "--roles", metaVar = "roles.in.subsystems", usage = "role definition file containing interesting roles")
    private File roleFile;

    @Argument(index = 0, metaVar = "inDir", usage = "input GTO directory")
    private File inDir;

    @Override
    protected void setDefaults() {
        this.roleFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        // Read the role map.
        if (this.roleFile == null)
            this.roleMap = new RoleMap();
        else {
            log.info("Reading interesting roles from {}.", this.roleFile);
            this.roleMap = RoleMap.load(this.roleFile);
        }
        // Create the maps.
        this.proteinMap = new HashMap<String, Set<Feature>>(2000000);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Create an MD5 computer.
        MD5Hex md5engine = new MD5Hex();
        // Get the input genomes.
        log.info("Scanning input directory {}.", this.inDir);
        GenomeDirectory genomes = new GenomeDirectory(this.inDir);
        log.info("{} genomes found in directory.", genomes.size());
        int gCount = 0;
        int fCount = 0;
        for (Genome genome : genomes) {
            log.info("Scanning {}.", genome);
            for (Feature feat : genome.getPegs()) {
                String seq = feat.getProteinTranslation();
                if (seq != null && ! seq.isEmpty()) {
                    Set<Feature> fSet = this.proteinMap.computeIfAbsent(md5engine.sequenceMD5(seq), x -> new HashSet<Feature>(10));
                    fSet.add(feat);
                    fCount++;
                }
            }
            gCount++;
            log.info("{} distinct proteins after {} genomes and {} features.", this.proteinMap.size(), gCount, fCount);
        }
        // Start the report.
        System.out.println("num\tfid\tfunction\tinteresting");
        // Now, for each sequence, we determine the most common function.
        FunctionMap funMap = new FunctionMap();
        CountMap<String> funCounts = new CountMap<String>();
        int protCount = 0;
        int badCount = 0;
        log.info("Producing report.");
        for (Set<Feature> flist : this.proteinMap.values()) {
            // Only bother if there is more than one feature with this protein.
            if (flist.size() > 1) {
                protCount++;
                funCounts.deleteAll();
                for (Feature feat : flist) {
                    Function fun = funMap.findOrInsert(feat.getPegFunction());
                    funCounts.count(fun.getId());
                }
                if (funCounts.size() > 1) {
                    // Here we have inconsistent annotation.
                    badCount++;
                    for (Feature feat : flist) {
                        boolean interest = feat.isInteresting(this.roleMap);
                        String iFlag = (interest ? "*" : "");
                        System.out.format("%8d\t%s\t%s\t%s%n", badCount, feat.getId(), feat.getPegFunction(), iFlag);
                    }
                    System.out.println();
                }
            }
        }
        log.info("{} inconsistent proteins found.  {} proteins occurred multiple times.", badCount, protCount);
    }

}
