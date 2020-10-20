/**
 *
 */
package org.theseed.genomes.compare;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.sequence.ProteinKmers;
import org.theseed.utils.BaseProcessor;

/**
 * This command uses kmers to compare identically-annotated features between close genomes.  If the features are
 * close enough, the gene name from the source feature is copied to the target feature.  The positional parameters
 * are the file names of the two genome IDs and the file name for the output genome.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	show more detailed progress messages
 * -m	maximum acceptable kmer distance for a match to be valid
 * -K	protein kmer size for distance computation
 *
 * @author Bruce Parrello
 *
 */
public class GeneCopyProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GeneCopyProcessor.class);
    /** function definition map */
    private FunctionMap funMap;
    /** source genome */
    private Genome source;
    /** target genome */
    private Genome target;
    /** map of function IDs to features */
    private Map<String, List<Feature>> funFeatures;
    /** map of feature IDs to aliases */
    private Map<String, String> aliasMap;

    // COMMAND-LINE OPTIONS

    /** maximum permissible distance */
    @Option(name ="-m", aliases = { "--maxDist" }, metaVar = "0.2", usage = "maximum permissible distance for a name transfer")
    private double maxDist;

    /** kmer size for distance measurement */
    @Option(name = "-K", aliases = { "--kmer", "--kmerSize" }, metaVar = "10", usage = "protein kmer size for distance computation")
    private int kmerSize;

    /** source genome */
    @Argument(index = 0, metaVar = "source.gto", usage = "source genome file", required = true)
    private File sourceFile;

    /** target genome */
    @Argument(index = 1, metaVar = "target.gto", usage = "genome file to update", required = true)
    private File targetFile;

    /** output file */
    @Argument(index = 2, metaVar = "output.gto", usage = "output file for modified genome", required = true)
    private File outputFile;


    @Override
    protected void setDefaults() {
        this.maxDist = 0.5;
        this.kmerSize = 8;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Validate the distance.
        if (this.maxDist < 0.0 || this.maxDist > 1.0)
            throw new IllegalArgumentException("Distance must be between 0 and 1.");
        // Validate the kmer size.
        if (this.kmerSize < 2)
            throw new IllegalArgumentException("Kmer size must be at least 2.");
        ProteinKmers.setKmerSize(this.kmerSize);
        // Validate the genomes.
        if (! this.sourceFile.canRead())
            throw new FileNotFoundException("Input genome file " + this.sourceFile + " not found or unreadable.");
        this.source = new Genome(this.sourceFile);
        if (! this.sourceFile.canRead())
            throw new FileNotFoundException("Input genome file " + this.targetFile + " not found or unreadable.");
        this.target = new Genome(this.targetFile);
        // Create the maps.
        this.funMap = new FunctionMap();
        this.funFeatures = new HashMap<String, List<Feature>>(4000);
        this.aliasMap = new HashMap<String, String>(4000);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Loop through the source, collecting features by role.
        log.info("Processing features in {}.", this.source);
        for (Feature feat : this.source.getPegs()) {
            Optional<String> geneName = feat.getAliases().stream()
                    .filter(x -> (! x.contains("|") && (x.length() == 3 || x.length() == 4))).findAny();
            if (geneName.isPresent()) {
                // Associate the feature with its function.
                String funDesc = feat.getPegFunction();
                Function fun = this.funMap.findOrInsert(funDesc);
                List<Feature> feats = this.funFeatures.computeIfAbsent(fun.getId(), x -> new ArrayList<Feature>(3));
                feats.add(feat);
                // Save the alias.
                this.aliasMap.put(feat.getId(), geneName.get());
            }
        }
        log.info("{} features with names, {} functions found.", this.aliasMap.size(), this.funFeatures.size());
        // Now loop through the target.
        int updates = 0;
        for (Feature feat : this.target.getPegs()) {
            // Find the best matching feature.
            String funDesc = feat.getPegFunction();
            Function fun = this.funMap.getByName(funDesc);
            if (fun != null) {
                List<Feature> feats = this.funFeatures.get(fun.getId());
                if (feats != null) {
                    // Here we have some features to check.  Find the closest.
                    ProteinKmers kmers = new ProteinKmers(feat.getProteinTranslation());
                    Feature found = null;
                    double fDist = this.maxDist;
                    for (Feature f2 : feats) {
                        ProteinKmers f2Kmers = new ProteinKmers(f2.getProteinTranslation());
                        double f2Dist = kmers.distance(f2Kmers);
                        if (f2Dist <= fDist) {
                            fDist = f2Dist;
                            found = f2;
                        }
                    }
                    if (found != null) {
                        String alias = this.aliasMap.get(found.getId());
                        log.debug("Feature {} with function \"{}\" has name {}.", feat.getId(), feat.getPegFunction(), alias);
                        feat.addAlias(alias);
                        updates++;
                    }
                }
            }
        }
        // Write out the target.
        log.info("Writing genome with {} updates to {}.", updates, this.outputFile);
        this.target.update(this.outputFile);
    }

}
