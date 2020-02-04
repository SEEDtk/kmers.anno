package org.theseed.proteins.kmers.anno;

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;

import org.apache.commons.codec.binary.StringUtils;
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

/**
 * This is the base class for annotating a genome.  It provides methods for handling the common
 * command-line processing and for annotating a single genome in place.
 *
 * @author Bruce Parrello
 *
 */
public class KmerProcessor {

    // FIELDS

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeKmerProcessor.class);

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

    // COMMAND LINE OPTIONS

    /** help option */
    @Option(name = "-h", aliases = { "--help" }, help = true)
    protected boolean help;
    /** minimum acceptable proposal strength */
    @Option(name = "-m", aliases = { "--minStrength",
            "--min" }, metaVar = "0.50", usage = "minimum acceptable proposal strength (0 to 1)")
    private double minStrength;
    /** maximum length factor */
    @Option(name = "-f", aliases = { "--fuzz", "--maxLength",
            "--max" }, metaVar = "2.0", usage = "maximum length increase factor for proteins (>= 1)")
    private double maxFuzz;
    /** contig kmer algorithm */
    @Option(name = "--algorithm", usage = "algorithm for retrieving contig kmers")
    private KmerFactory.Type kmerType;

    /** kmer length */
    @Option(name = "-K", aliases = { "--kmer" }, metaVar = "10", usage = "protein kmer length (default 6)")
    private void setKmers(int newSize) {
        KmerReference.setKmerSize(newSize);
    }

    /** number of close genomes to use */
    @Option(name = "-n", aliases = { "--nGenomes",
            "--num" }, metaVar = "2", usage = "maximum number of close genomes to scan")
    private int maxGenomes;
    @Option(name = "--trace", usage = "function assignment to be traced")
    private String traceFunction;

    public KmerProcessor() {
        super();
    }

    /**
     * Verify that the command-line options are correct.
     */
    protected void validateParms() {
        // Verify the options.
        if (this.minStrength >= 1.0)
            throw new IllegalArgumentException("Minimum strength must be less than 1.");
        if (this.maxFuzz <= 1.0)
            throw new IllegalArgumentException("Length fuzz factor must be greater than 1.");
        // Create the kmer factory.
        this.kmerFactory = KmerFactory.create(kmerType);
        // Connect to PATRIC.
        log.info("Connecting to PATRIC.");
        this.p3 = new Connection();
    }

    /**
     * Set the default command-line options.
     */
    protected void setDefaults() {
        this.minStrength = 0.20;
        this.maxFuzz = 1.5;
        this.maxGenomes = 10;
        this.kmerType = KmerFactory.Type.AGGRESSIVE;
    }

    /**
     * Annotate a genome.
     *
     * @param genome		genome to annotate
     */
    protected void annotateGenome(Genome genome) {
        log.info("Annotating proposed genome {}: {}", genome.getId(), genome.getName());
        // Our proposed annotations will be accumulated in here.
        log.info("Minimum proposal strength is {}.  Kmer size is {}.", this.minStrength, KmerReference.getKmerSize());
        // Note that because the evidence is protein evidence, and the proposal is DNA, we need to
        // divide the strength by 3 to make it work.  The user will never see the real strength.
        double realStrength = this.minStrength / 3;
        PegProposalList proposals = new PegProposalList(genome, realStrength);
        // Get the contig kmers from the new genome.
        Map<String, Collection<Location>> contigKmers = this.kmerFactory.findKmers(genome);
        log.info("{} kmers found in genome.", contigKmers.size());
        // Extract the close genomes.
        SortedSet<CloseGenome> closeGenomes = genome.getCloseGenomes();
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
                            PegProposal found = proposals.propose(wholeLoc, peg.getFunction(), evidenceCount);
                            // Here we check for the function we're debugging.  This equals function is null-safe.
                            if (found != null && StringUtils.equals(this.traceFunction, peg.getFunction())) {
                                log.info("Proposal stored using {} at location {} with evidence {} and strength {}.", peg.getId(), wholeLoc,
                                        evidenceCount, found.getStrength());
                            }
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
        // Now we want to loop through the proposals, creating features.
        // Initialize the peg counter.  We use this to generate IDs.
        this.pegCount = 0;
        // Get a DNA translator.
        this.xlator = new DnaTranslator(genome.getGeneticCode());
        log.info("Using genetic code {}.", genome.getGeneticCode());
        for (PegProposal current : proposals) {
            this.makeFeature(current, genome);
        }
        // All done.
        log.info("Processing complete. {} features in genome.", this.pegCount);
    }

    /**
     * Store a peg proposal as a feature in the genome.
     *
     * @param proposal	location and function of the proposed peg
     * @param genome	genome to contain the new feature
     */
    private void makeFeature(PegProposal proposal, Genome genome) {
        // Create the feature ID.
        this.pegCount++;
        String fid = String.format("fig|%s.peg.%d", genome.getId(), this.pegCount);
        // Get the proposed location.
        Location loc = proposal.getLoc();
        // Create the feature.
        Feature feat = new Feature(fid, proposal.getFunction(), loc.getContigId(), loc.getStrand(), loc.getLeft(), loc.getRight());
        // Compute the protein translation.  Note we have to trim the stop codon off the translation.
        String dna = genome.getDna(loc);
        String prot = this.xlator.pegTranslate(dna, 1, dna.length() - 3);
        feat.setProteinTranslation(prot);
        // Store the feature in the genome.
        genome.addFeature(feat);
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
