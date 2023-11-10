/**
 *
 */
package org.theseed.genome.compare;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.sequence.ProteinKmers;

/**
 * This command uses kmers to compare identically-annotated features between close genomes.  If the features are
 * close enough, the aliases from the source feature are copied to the target feature.  The positional parameters
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
    private Map<String, Collection<String>> aliasMap;

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
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Validate the distance.
        if (this.maxDist < 0.0 || this.maxDist > 1.0)
            throw new ParseFailureException("Distance must be between 0 and 1.");
        // Validate the kmer size.
        if (this.kmerSize < 2)
            throw new ParseFailureException("Kmer size must be at least 2.");
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
        this.aliasMap = new HashMap<String, Collection<String>>(4000);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Loop through the source, collecting features by role.
        log.info("Processing features in {}.", this.source);
        for (Feature feat : this.source.getPegs()) {
            Collection<String> aliases = feat.getAliases();
            if (aliases != null && aliases.size() > 0) {
                // Associate the feature with its function.
                String funDesc = feat.getPegFunction();
                Function fun = this.funMap.findOrInsert(funDesc);
                List<Feature> feats = this.funFeatures.computeIfAbsent(fun.getId(), x -> new ArrayList<Feature>(3));
                feats.add(feat);
                // Save the aliases.
                this.aliasMap.put(feat.getId(), aliases);
            }
        }
        log.info("{} features with aliases, {} functions found.", this.aliasMap.size(), this.funFeatures.size());
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
                        Collection<String> aliases = this.aliasMap.get(found.getId());
                        if (log.isDebugEnabled())
                            log.debug("Feature {} with function \"{}\" has aliases {}.", feat.getId(), feat.getPegFunction(),
                                    aliases.stream().collect(Collectors.joining(", ")));
                        aliases.stream().forEach(x -> feat.addAlias(x));
                        updates++;
                    }
                }
            }
        }
        // Write out the target.
        log.info("Writing genome with {} updates to {}.", updates, this.outputFile);
        this.target.save(this.outputFile);
    }

}
