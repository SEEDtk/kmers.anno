/**
 *
 */
package org.theseed.locations;

import org.theseed.genome.Genome;

/**
 * This class represents a proposed protein annotation in an ORF.  Each proposal can be put
 * into a TreeMap or HashMap.
 *
 * @author Bruce Parrello
 *
 */
public class PegProposal implements Comparable<PegProposal> {

    // FIELDS

    /** location of the protein (includes the stop) */
    private Location loc;
    /** proposed functional assignment */
    private String function;
    /** evidence for this proposal */
    private double strength;


    /**
     * Create a new peg proposal.
     *
     * @param loc		location containing the peg
     * @param function	functional assignment for this peg
     * @param strength	evidence for this proposal
     */
    private PegProposal(Location loc, String function, double strength) {
        this.loc = loc;
        this.function = function;
        this.strength = strength;
    }

    /**
     * Attempt to propose a peg for a genome.
     *
     * @param genome	genome of interest
     * @param loc		location in the genome of the proposed peg
     * @param function	functional assignment for this peg
     * @param evidence	number of base pairs of evidence for this proposal
     *
     * @return the peg proposal, fully extended to a start and a stop, or NULL if the proposal is invalid
     */
    public static PegProposal create(Genome genome, Location loc, String function, int evidence) {
        PegProposal retVal = null;
        // Try extending the location.
        Location realLocation = loc.extend(genome);
        if (realLocation != null) {
            double strength = ((double) evidence) / realLocation.getLength();
            retVal = new PegProposal(realLocation, function, strength);
        }
        return retVal;
    }

    /**
     * Override this peg proposal with the information in another.  Note that if the locations
     * are not compatible, the results will be ridiculous.
     *
     * @param other		new proposal to merge
     */
    public void merge(PegProposal other) {
        this.function = other.function;
        this.loc.setBegin(other.loc.getBegin());
        this.strength = other.strength;
    }



    /**
     * Compare two peg proposals for sorting.  They are equal if they have the same end point and strand.
     *
     * Otherwise, they are ordered inside the contig by left edge and then length.
     *
     * @param other peg proposal to compare
     *
     * @return a negative value if we should sort earlier, a positive value if we should sort later,
     * 		   0 if we are equal
     */
    @Override
    public int compareTo(PegProposal other) {
        // Sort by contig first.
        int retVal = this.loc.getContigId().compareTo(other.loc.getContigId());
        if (retVal == 0) {
            // Only proceed if the endpoints or strands are different.
            if (this.loc.getEnd() != other.loc.getEnd() || this.loc.getDir() != other.loc.getDir()) {
                retVal = this.loc.getLeft() - other.loc.getLeft();
                if (retVal == 0) {
                    // Shorter locations go first.
                    retVal = this.loc.getLength() - other.loc.getLength();
                }
            }
        }
        return retVal;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        if (loc != null) {
            result = result + loc.getContigId().hashCode();
            result = result * prime + loc.getDir();
            result = result * prime + loc.getEnd();
        }
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof PegProposal)) {
            return false;
        }
        Location other = ((PegProposal) obj).loc;
        return (other.getContigId() == this.loc.getContigId() && other.getEnd() == this.loc.getEnd() && other.getDir() == this.loc.getDir());
    }

    /**
     * @return the proposed peg location
     */
    public Location getLoc() {
        return loc;
    }

    /**
     * @return the proposed functional assignment
     */
    public String getFunction() {
        return function;
    }

    /**
     * @return the strength of this proposal's evidence
     */
    public double getStrength() {
        return strength;
    }

    /**
     * @return TRUE if this proposal is superior to a specified other proposal.
     *
     * @param oldProposal	other proposal for comparison
     */
    public boolean betterThan(PegProposal oldProposal) {
        boolean retVal = (this.strength > oldProposal.strength);
        if (! retVal && this.strength == oldProposal.strength)
            retVal = (this.getLoc().getLength() > oldProposal.getLoc().getLength());
        return retVal;
    }

}
