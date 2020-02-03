/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.genome.CloseGenome;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.FramedLocationLists;
import org.theseed.locations.Location;
import org.theseed.locations.PegProposal;
import org.theseed.locations.PegProposalList;
import org.theseed.locations.SortedLocationList;
import org.theseed.p3api.Connection;
import org.theseed.p3api.P3Genome;
import org.theseed.p3api.P3Genome.Details;
import org.theseed.proteins.DnaTranslator;
import org.theseed.proteins.kmers.KmerFactory;
import org.theseed.proteins.kmers.KmerReference;
import org.theseed.utils.ICommand;

/**
 * This program attempts to project proteins from close genomes to a new genome.  The new genome should have its taxonomy
 * information, close-genome list, genetic code, ID, and name filled in, as well as all of its contigs.  We loop through
 * the N closest genomes, using kmers to project proteins onto ORFs.
 *
 * The genome should be in the form of a CTO.  It is read from the standard input, and an annotated GTO is written to the
 * standard output.
 *
 * The following command-line options are supported.
 *
 * -m	minimum strength for an annotation proposal to be acceptable
 * -f	maximum factor for length increase in forming a proposal
 * -K	kmer length to use
 * -n	number of close genomes to use
 * -i	name of the input file (overrides STDIN)
 * -o	name of the output file (overrides STDOUT)
 *
 * --compare	if specified, the name of a GTO file for a genome with the same contigs; the new genome will be compared to
 * 				the genome in the GTO file
 *
 * @author Bruce Parrello
 *
 */
public class KmerProcessor implements ICommand {

    // FIELDS

    /** input genome */
    private Genome newGenome;
    /** connection to PATRIC */
    private Connection p3;
    /** list of locations for each feature, separated vby frame */
    private FramedLocationLists framer;
    /** number of pegs created */
    private int pegCount;
    /** DNA translator */
    private DnaTranslator xlator;
    /** contig kmer factory */
    private KmerFactory kmerFactory;

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(KmerProcessor.class);

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** minimum acceptable proposal strength */
    @Option(name = "-m", aliases = { "--minStrength", "--min" }, metaVar = "0.50", usage = "minimum acceptable proposal strength (0 to 1)")
    private double minStrength;

    /** maximum length factor */
    @Option(name = "-f", aliases = { "--fuzz", "--maxLength", "--max" }, metaVar = "2.0", usage = "maximum length increase factor for proteins (>= 1)")
    private double maxFuzz;

    /** kmer length */
    @Option(name = "-K", aliases = { "--kmer" }, metaVar = "10", usage = "protein kmer length")
    private void setKmers(int newSize) {
        KmerReference.setKmerSize(newSize);
    }

    /** contig kmer algorithm */
    @Option(name = "--algorithm", usage = "algorithm for retrieving contig kmers")
    private KmerFactory.Type kmerType;

    /** number of close genomes to use */
    @Option(name = "-n", aliases = { "--nGenomes", "--num" }, metaVar = "2", usage = "maximum number of close genomes to scan")
    private int maxGenomes;

    /** input file name */
    @Option(name = "-i", aliases = { "--input" }, usage = "input file name (if not STDIN)")
    private File inFile;

    /** output file name */
    @Option(name = "-o", aliases = { "--output" }, usage = "output file name (if not STDOUT)")
    private File outFile;


    /**
     * Parse command-line options to specify the parameters of this object.
     *
     * @param args	an array of the command-line parameters and options
     *
     * @return TRUE if successful, FALSE if the parameters are invalid
     */
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.help = false;
        this.minStrength = 0.20;
        this.maxFuzz = 1.5;
        this.maxGenomes = 10;
        this.inFile = null;
        this.outFile = null;
        this.kmerType = KmerFactory.Type.AGGRESSIVE;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                // Verify the options.
                if (this.minStrength >= 1.0)
                    throw new IllegalArgumentException("Minimum strength must be less than 1.");
                if (this.maxFuzz <= 1.0)
                    throw new IllegalArgumentException("Length fuzz factor must be greater than 1.");
                // Create the kmer factory.
                this.kmerFactory = KmerFactory.create(kmerType);
                // Denote we can run.
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            // For parameter errors, we display the command usage.
            parser.printUsage(System.err);
        }
        return retVal;
    }

    public void run() {
        try {
            log.info("Connecting to PATRIC.");
            this.p3 = new Connection();
            if (this.inFile != null) {
                log.info("Reading genome from {}.", this.inFile);
                this.newGenome = new Genome(this.inFile);
            } else {
                log.info("Reading genome from standard input.");
                this.newGenome = new Genome(System.in);
            }
            log.info("Annotating proposed genome {}: {}", this.newGenome.getId(), this.newGenome.getName());
            // Our proposed annotations will be accumulated in here.
            log.info("Minimum proposal strength is {}.", this.minStrength);
            // Note that because the evidence is protein evidence, and the proposal is DNA, we need to
            // divide the strength by 3 to make it work.  The user will never see the real strength.
            double realStrength = this.minStrength / 3;
            PegProposalList proposals = new PegProposalList(this.newGenome, realStrength);
            // Get the contig kmers from the new genome.
            Map<String, Collection<Location>> contigKmers = this.kmerFactory.findKmers(this.newGenome);
            log.info("{} kmers found in genome.", contigKmers.size());
            // Extract the close genomes.
            SortedSet<CloseGenome> closeGenomes = this.newGenome.getCloseGenomes();
            log.info("{} close genomes available from input.", closeGenomes.size());
            // Create the framed location lists.
            this.framer = new FramedLocationLists();
            // Iterate through the close genomes, stopping at the maximum.
            Iterator<CloseGenome> iter = closeGenomes.iterator();
            int iGenome = 1;
            while (iter.hasNext() && iGenome <= this.maxGenomes) {
                CloseGenome closeGenome = iter.next();
                log.info("Retrieving close genome #{} {}: {}.", iGenome,
                        closeGenome.getGenomeId(), closeGenome.getGenomeName());
                Genome oldGenome = P3Genome.Load(p3, closeGenome.getGenomeId(), Details.PROTEINS);
                if (oldGenome == null) {
                    log.warn("Genome {} not found-- skipping.");
                } else {
                    iGenome++;
                    // Get the peg kmers for this genome.
                    Set<KmerReference> pegKmers = this.getPegKmers(oldGenome);
                    // Match the two kmer sets to build the location lists.
                    for (KmerReference pegKmer : pegKmers) {
                        // Find where this kmer occurs in the new genome.
                        Collection<Location> kmerLocs = contigKmers.get(pegKmer.getKmer());
                        if (kmerLocs != null) {
                            // We have a match.  Add these as possible locations for the associated peg.
                            // Note that the "contigID()" in the peg kmer is actually the peg ID, so
                            // the contig location (kmerLoc) will be connected to a peg ID.
                            for (Location kmerLoc : kmerLocs)
                                this.framer.connect(pegKmer.getLoc().getContigId(), kmerLoc);
                        }
                    }
                    log.info("{} matching kmers found.", this.framer.size());
                    // Now we loop through the location lists, proposing pegs.
                    int pegsFound = 0;
                    int lowKmerCountFound = 0;
                    int proposalCount = 0;
                    for (FramedLocationLists.Report report : this.framer) {
                        SortedLocationList locList = report.getList();
                        // Get the feature corresponding to this report.
                        Feature peg = oldGenome.getFeature(report.getId());
                        pegsFound++;
                        // Determine the minimum number of kmers required to have sufficient
                        // strength and the maximum proposal length.  Note that we scale
                        // the protein length to base pairs, because all our locations
                        // are DNA locations.
                        int pegLen = peg.getProteinLength() * 3;
                        int maxLen = (int) (pegLen * this.maxFuzz + 1);
                        int minKmers = (int) (pegLen * realStrength);
                        // Insure we have enough matches to even look.
                        if (minKmers > locList.size())
                            lowKmerCountFound++;
                        else {
                            // Determine the last location that can possibly start a proposal and
                            // loop through all the possible starts.
                            int n = locList.size() - minKmers;
                            for (int i = 0; i <= n; i++) {
                                Location firstLoc = locList.get(i);
                                // The evidence count starts at 1 because it includes the first location.
                                int evidenceCount = 1;
                                // Determine the last permissible right edge and the actual right
                                // edge.
                                int maxEdge = firstLoc.getLeft() + maxLen;
                                int bestEdge = firstLoc.getRight();
                                for (Location loc : locList.contigRange(i)) {
                                    if (loc.getRight() < maxEdge) {
                                        evidenceCount++;
                                        bestEdge = Math.max(loc.getRight(), bestEdge);
                                    }
                                }
                                // Propose this peg.
                                Location wholeLoc = Location.create(firstLoc.getContigId(), firstLoc.getStrand(), firstLoc.getLeft(), bestEdge);
                                proposals.propose(wholeLoc, peg.getFunction(), evidenceCount);
                                proposalCount++;
                            }
                        }
                    }
                    log.info("{} peg/frame pairs examined, {} had too few kmers, {} proposals were made.",
                            pegsFound, lowKmerCountFound, proposalCount);
                }
                // Clear the location lists for the next pass.
                framer.clear();
            }
            // Write the proposal statistics.  Note that we do not list the proposals
            // discarded because they were weaker versions of existing proposals.
            log.info("{} proposals made, {} merged, {} rejected, {} too weak, {} kept.",
                    proposals.getMadeCount(), proposals.getMergeCount(), proposals.getRejectedCount(),
                    proposals.getWeakCount(), proposals.getProposalCount());
            // Now we want to loop through the proposals. The proposals are presented in location order.  We keep one proposal in
            // reserve.  If the next proposal overlaps, we keep the strongest.  If it does not, we convert the proposal
            // to a feature in the new genome.
            int overlapCount = 0;
            Iterator<PegProposal> pIter = proposals.iterator();
            // Error out if there are no proposals.
            if (! pIter.hasNext())
                throw new RuntimeException("No matching proteins found.  Unable to annotate this genome.");
            // Initialize the peg counter.  We use this to generate IDs.
            this.pegCount = 0;
            // Get a DNA translator.
            this.xlator = new DnaTranslator(newGenome.getGeneticCode());
            log.info("Using genetic code {}.", newGenome.getGeneticCode());
            // Loop through the proposals.
            PegProposal reserve = pIter.next();
            Location reserveLoc = reserve.getLoc();
            while (pIter.hasNext()) {
                PegProposal current = pIter.next();
                if (reserveLoc.distance(current.getLoc()) < 0) {
                    // Here the proposals overlap.  This generally means they are on different frames.  We keep the
                    // better one.
                    overlapCount++;
                    if (current.betterThan(reserve)) {
                        reserve = current;
                        reserveLoc = current.getLoc();
                    }
                } else {
                    // Here the proposals do not overlap.  Output the reserve, and save the current one.
                    this.makeFeature(reserve);
                    reserve = current;
                    reserveLoc = current.getLoc();
                }
            }
            // All done.  Make a feature from the last reserve.
            this.makeFeature(reserve);
            // Write the result.
            if (this.outFile != null) {
                log.info("Writing genome to {}.", this.outFile);
                this.newGenome.update(outFile);
            } else {
                log.info("Writing genome to standard output.");
                this.newGenome.update(System.out);
            }
            log.info("Processing complete. {} features in genome. {} overlaps discarded", this.pegCount, overlapCount);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Store a peg proposal as a feature in the genome.
     *
     * @param proposal	location and function of the proposed peg.
     */
    private void makeFeature(PegProposal proposal) {
        // Create the feature ID.
        this.pegCount++;
        String fid = String.format("fig|%s.peg.%d", this.newGenome.getId(), this.pegCount);
        // Get the proposed location.
        Location loc = proposal.getLoc();
        // Create the feature.
        Feature feat = new Feature(fid, proposal.getFunction(), loc.getContigId(), loc.getStrand(), loc.getLeft(), loc.getRight());
        // Compute the protein translation.  Note we have to trim the stop codon off the translation.
        String dna = this.newGenome.getDna(loc);
        String prot = this.xlator.pegTranslate(dna, 1, dna.length() - 3);
        feat.setProteinTranslation(prot);
        // Store the feature in the genome.
        this.newGenome.addFeature(feat);
    }

    /**
     * @return a set of the unique kmers in the specified genome's pegs.
     *
     * @param genome		genome of interest
     */
    private Set<KmerReference> getPegKmers(Genome genome) {
        // Get all the kmers.
        log.info("Scanning for kmers in {} {}.", genome.getId(), genome.getName());
        CountMap<KmerReference> pegKmers = KmerReference.countPegKmers(genome);
        // Keep the unique ones.
        Set<KmerReference> retVal = pegKmers.getSingletons();
        log.info("{} total protein kmers found in the features. {} were unique.", pegKmers.size(), retVal.size());
        return retVal;
    }

}
