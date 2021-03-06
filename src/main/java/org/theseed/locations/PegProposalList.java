/**
 *
 */
package org.theseed.locations;

import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.theseed.genome.Genome;

/**
 * Manage a list of peg proposals for a genome.  This object will only accept proposals with a specified
 * strength, and only keep the strongest proposal for a given ORF.  It is iterable in contig order.  That
 * is, the natural order of this object will present the proposed pegs in the way they should be numbered.
 *
 * @author Bruce Parrello
 *
 */
public class PegProposalList implements Iterable<PegProposal> {

    // FIELDS
    /** minimum acceptable strength */
    private double minStrength;
    /** number of proposals made */
    private int madeCount;
    /** number of invalid proposals */
    private int rejectedCount;
    /** number of weak proposals */
    private int weakCount;
    /** number of insufficient-evidence proposals */
    private int smallCount;
    /** number of proposals merged */
    private int mergeCount;
    /** list of proposals */
    SortedSet<PegProposal> proposals;
    /** target genome for the proposals */
    Genome genome;
    /** minimum acceptable evidence */
    private int minEvidence;

    /**
     * Create a new peg proposal list.
     */
    public PegProposalList(Genome genome, double minStrength, int minEvidence) {
        this.genome = genome;
        this.minStrength = minStrength;
        this.minEvidence = minEvidence;
        this.madeCount = 0;
        this.rejectedCount = 0;
        this.mergeCount = 0;
        this.weakCount = 0;
        this.smallCount = 0;
        this.proposals = new TreeSet<PegProposal>();
    }

    /**
     * Propose a peg.
     *
     * @param loc		the location for the proposal
     * @param function	the proposed functional assignment
     * @param evidence	the number of base pairs of evidence for the proposal; the strength is
     * 					the evidence divided by the final location length
     *
     * @return the new proposal if it took, else FALSE
     */
    public PegProposal propose(Location loc, String function, int evidence) {
        PegProposal retVal = null;
        this.madeCount++;
        PegProposal newProposal = PegProposal.create(this.genome, loc, function, evidence);
        if (newProposal == null) {
            this.rejectedCount++;
        } else if (newProposal.getStrength() < this.minStrength) {
            this.weakCount++;
        } else if (evidence < this.minEvidence) {
            this.smallCount++;
        } else {
            // The proposal is worth keeping.  Try to add it. TRUE means it was added.
            boolean added = this.proposals.add(newProposal);
            if (added) {
                retVal = newProposal;
            } else {
                // Here there is a duplicate proposal.  Compare it with the current one.
                PegProposal oldProposal = this.proposals.tailSet(newProposal).first();
                if (newProposal.betterThan(oldProposal)) {
                    oldProposal.merge(newProposal);
                    this.mergeCount++;
                    retVal = oldProposal;
                }
            }
        }
        return retVal;
    }

    /**
     * @return the number of proposals made
     */
    public int getMadeCount() {
        return madeCount;
    }

    /**
     * @return the number of proposals rejected for insufficient strength
     */
    public int getRejectedCount() {
        return rejectedCount;
    }

    /**
     * @return the number of proposals merged into weaker ones
     */
    public int getMergeCount() {
        return mergeCount;
    }

    /**
     * @return the number of active proposals
     */
    public int getProposalCount() {
        return this.proposals.size();
    }

    @Override
    public Iterator<PegProposal> iterator() {
        return this.proposals.iterator();
    }

    /**
     * @return the number of proposals discarded for being too weak
     */
    public int getWeakCount() {
        return weakCount;
    }

    /**
     * @return the number of proposals discarded for too little evidence
     */
    public int getSmallCount() {
        return smallCount;
    }

}
