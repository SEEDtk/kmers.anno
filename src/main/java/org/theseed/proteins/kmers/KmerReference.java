/**
 *
 */
package org.theseed.proteins.kmers;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.HashMap;

import org.apache.commons.lang3.StringUtils;
import org.theseed.counters.CountMap;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.DnaTranslator;

/**
 * This object represents a protein kmer at a specific location.  It is keyed on the kmer
 * text, so two identical kmers at different locations will compare equal.
 *
 * @author Bruce Parrello
 *
 */
public class KmerReference {

    /** protein kmer size */
    private static int K = 6;

    /** offset from start of kmer to end */
    private static int K3 = 17;

    // FIELDS

    /** kmer sequence */
    String kmer;
    /** location of the kmer */
    Location loc;

    /**
     * Change the kmer size.
     *
     * @param newSize	proposed new kmer size
     */
    public static void setKmerSize(int newSize) {
        K = newSize;
        K3 = K * 3 - 1;
    }

    /**
     * @return the kmer size
     */
    public static int getKmerSize() {
        return K;
    }

    /**
     * Create a reference to a specific kmer at a specific location.
     *
     * @param kmer		protein kmer sequence
     * @param seqId		ID of the sequence containing the kmer
     * @param begin		start position of the kmer (1-based)
     * @param dir		strand of the kmer (+ or -)
     */
    public KmerReference(String kmer, String seqId, int begin, String dir) {
        // Create the location.
        this.loc = Location.create(seqId, dir, begin, begin + K3);
        // Store the kmer string.
        this.kmer = kmer;
    }

    /**
     * @return the hash code for this object; note we only hash the kmer
     */
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((kmer == null) ? 0 : kmer.hashCode());
        return result;
    }

    /**
     * @return TRUE if this kmer reference is for the same kmer as another reference
     */
    @Override
    public boolean equals(Object obj) {
        boolean retVal;
        if (this == obj)
            retVal = true;
        else if (!(obj instanceof KmerReference))
            retVal = false;
        else {
            KmerReference other = (KmerReference) obj;
            if (kmer == null)
                retVal = (other.kmer == null);
            else
                retVal = kmer.equals(other.kmer);
        }
        return retVal;
    }

    /**
     * @return the kmer sequence
     */
    public String getKmer() {
        return kmer;
    }

    /**
     * @return the location of the kmer
     */
    public Location getLoc() {
        return loc;
    }

    /** Count all the kmers in the pegs of the specified genome.  Each kmer will be associated with its first occurrence.
     *
     * @param genome	genome whose peg kmers are to be counted
     *
     * @return a count map of the kmers in the genome's protein-encoding genes
     */
    public static CountMap<KmerReference> countPegKmers(Genome genome) {
        CountMap<KmerReference> retVal = new CountMap<KmerReference>();
        // Loop through the pegs.
        for (Feature feat : genome.getPegs()) {
            // Save the feature ID.
            String fid = feat.getId();
            // Get the protein translation.
            String prot = feat.getProteinTranslation();
            if (prot != null) {
                // Compute the location of the last kmer.
                int end = prot.length() - K;
                // Loop through the protein.
                for (int i = 0; i < end; i++) {
                    String kmer = prot.substring(i, i + K);
                    // Note we make sure there are no ambiguity characters.
                    if (StringUtils.containsNone(kmer, 'X')) {
                        KmerReference kmerRef = new KmerReference(kmer, fid, i + 1, "+");
                        retVal.count(kmerRef);
                    }
                }
            }
        }
        return retVal;
    }

    /**
     * Get all the kmers in the contigs of the specified genome.  Each kmer will be associated with its a list
     * of locations.
     *
     * @param genome	genome whose contig kmers are to be counted
     *
     * @return a map of the kmer locations in the genome
     */
    public static Map<String, Collection<Location>> getContigKmers(Genome genome) {
        Map<String, Collection<Location>> retVal = new HashMap<String, Collection<Location>>(genome.getLength() * 2);
        // Get a DNA translator.
        DnaTranslator xlator = new DnaTranslator(genome.getGeneticCode());
        // Loop through the contigs.
        for (Contig contig : genome.getContigs()) {
            // Count the kmers on the forward strand.
            processKmers(retVal, xlator, contig.getId(), new KmerPosition.Plus(contig), contig.getSequence());
            // Count the kmers on the reverse strand.
            processKmers(retVal, xlator, contig.getId(), new KmerPosition.Minus(contig), contig.getRSequence());
        }
        return retVal;
    }

    /**
     * Count the kmers in a single strand of a contig.
     *
     * @param counters		kmer counting map
     * @param xlator		DNA translator for converting to proteins
     * @param contigId		ID of the source contig
     * @param calculator	position calculator for this strand
     * @param sequence		DNA sequence of the strand
     */
    private static void processKmers(Map<String, Collection<Location>> kmers, DnaTranslator xlator, String contigId,
            KmerPosition calculator, String sequence) {
        // Here we go through each frame, extracting kmers.  We translate the entire sequence, then pull substrings.
        for (int frame = 1; frame <= 3; frame++) {
            String frameProtein = xlator.translate(sequence, frame, sequence.length());
            // Compute the point past the last legal kmer start.
            int end = frameProtein.length() - K;
            for (int i = 0; i < end; i++) {
                // Get the kmer.  We ignore it if it contains ambiguity characters or a stop.
                String kmerString = frameProtein.substring(i, i + K);
                if (StringUtils.containsNone(kmerString, '*', 'X')) {
                    int left = calculator.calcLeft(i, frame);
                    Location loc = Location.create(contigId, calculator.getDir(), left, left + K3);
                    // Associate the location with the kmer.
                    Collection<Location> locs = kmers.get(kmerString);
                    if (locs == null) {
                        locs = new ArrayList<Location>(5);
                        kmers.put(kmerString, locs);
                    }
                    locs.add(loc);
                }
            }
        }
    }

    @Override
    public String toString() {
        return "KmerReference [" + (kmer != null ? "kmer=" + kmer + ", " : "") + (loc != null ? "loc=" + loc : "")
                + "]";
    }

}
