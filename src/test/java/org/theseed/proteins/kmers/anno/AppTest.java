package org.theseed.proteins.kmers.anno;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Frame;
import org.theseed.locations.FramedLocationLists;
import org.theseed.locations.Location;
import org.theseed.locations.PegProposal;
import org.theseed.locations.PegProposalList;
import org.theseed.locations.SortedLocationList;
import org.theseed.proteins.CodonSet;
import org.theseed.proteins.DnaTranslator;
import org.theseed.proteins.kmers.KmerFactory;
import org.theseed.proteins.kmers.KmerReference;

/**
 * Unit test for simple App.
 */
public class AppTest extends TestCase
{
    private static final String TEST_CONTIG = "51203.13.con.0001";

    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    // SAMPLE KMERS

    private static final String KMER1 = "SNENIQNA";
    private static final String KMER2 = "SNENMQNA";
    private static final String KMER3 = "GDQM";

    /**
     * basic kmer test
     */
    public void testKmerReferences() {
        KmerReference.setKmerSize(8);
        KmerReference k1 = new KmerReference(KMER1, "contig1", 10, "+");
        assertThat(k1.getKmer(), equalTo(KMER1));
        Location loc = k1.getLoc();
        assertThat(loc.getContigId(), equalTo("contig1"));
        assertThat(loc.getLeft(), equalTo(10));
        assertThat(loc.getLength(), equalTo(24));
        assertThat(loc.getBegin(), equalTo(10));
        assertThat(loc.getDir(), equalTo('+'));
        KmerReference k1a = new KmerReference(KMER1, "contig1", 100, "-");
        assertThat(k1a.getKmer(), equalTo(KMER1));
        assertTrue(k1.equals(k1a));
        assertThat(k1.hashCode(), equalTo(k1a.hashCode()));
        loc = k1a.getLoc();
        assertThat(loc.getContigId(), equalTo("contig1"));
        assertThat(loc.getLeft(), equalTo(100));
        assertThat(loc.getLength(), equalTo(24));
        assertThat(loc.getEnd(), equalTo(100));
        assertThat(loc.getDir(), equalTo('-'));
        KmerReference k2 = new KmerReference(KMER2, "contig1", 100, "-");
        assertThat(k2.getKmer(), equalTo(KMER2));
        assertFalse(k1.equals(k2));
        Location loc2 = k2.getLoc();
        assertThat(loc2, equalTo(loc));
        KmerReference.setKmerSize(4);
        KmerReference k3 = new KmerReference(KMER3, "contig2", 10, "+");
        assertThat(k3.getKmer(), equalTo(KMER3));
        loc = k3.getLoc();
        assertThat(loc.getContigId(), equalTo("contig2"));
        assertThat(loc.getLeft(), equalTo(10));
        assertThat(loc.getBegin(), equalTo(10));
        assertThat(loc.getLength(), equalTo(12));
    }

    /**
     * test kmer counting
     *
     * @throws IOException
     */
    public void testKmerContigCounts() throws IOException {
        KmerReference.setKmerSize(8);
        Genome smallGto = new Genome(new File("src/test", "small.gto"));
        DnaTranslator xlator = new DnaTranslator(smallGto.getGeneticCode());
        Map<String, Collection<Location>> kmerLocs = KmerReference.getContigKmers(smallGto);
        // Verify that we found the kmers in all the right places.
        verifyKmers(smallGto, xlator, kmerLocs);
        KmerFactory factory = KmerFactory.create(KmerFactory.Type.AGGRESSIVE);
        kmerLocs = factory.findKmers(smallGto);
        verifyKmers(smallGto, xlator, kmerLocs);
        factory = KmerFactory.create(KmerFactory.Type.STRICT);
        kmerLocs = factory.findKmers(smallGto);
        verifyKmers(smallGto, xlator, kmerLocs);
        for (String kmer : kmerLocs.keySet()) {
            Collection<Location> locs = kmerLocs.get(kmer);
            assertThat(locs.size(), equalTo(1));
        }
    }

    /**
     * Verify that a kmer/location map is valid
     */
    private void verifyKmers(Genome genome, DnaTranslator xlator, Map<String, Collection<Location>> kmerLocs) {
        for (String kmer : kmerLocs.keySet()) {
            for (Location loc : kmerLocs.get(kmer)) {
                assertFalse("Invalid characters in " + kmer, StringUtils.containsAny(kmer, 'X', '*'));
                assertThat("Kmer does not match location " + loc, kmer, equalTo(xlator.translate(genome.getDna(loc))));
            }
        }
    }

    /**
     * test kmer peg counting
     *
     * @throws IOException
     */
    public void testKmerPegCounts() throws IOException {
        KmerReference.setKmerSize(8);
        Genome smallGto = new Genome(new File("src/test", "small.gto"));
        CountMap<KmerReference> kmerCounts = KmerReference.countPegKmers(smallGto);
        // Verify that we found the kmers in all the right places.
        for (CountMap<KmerReference>.Count kmerCount : kmerCounts.counts()) {
            KmerReference kmerRef = kmerCount.getKey();
            String kmer = kmerRef.getKmer();
            Location loc = kmerRef.getLoc();
            assertFalse("Invalid characters in " + kmerRef, StringUtils.containsAny(kmer, 'X', '*'));
            // Get the peg's protein.
            Feature feat = smallGto.getFeature(loc.getContigId());
            String prot = feat.getProteinTranslation();
            int offset = loc.getLeft() - 1;
            assertThat("Kmer does not match location in " + kmerRef, kmer, equalTo(prot.substring(offset, offset + 8)));
        }
    }

    /**
     * test peg proposals
     *
     * @throws IOException
     */
    public void testPegProposals() throws IOException {
        CodonSet starts = new CodonSet("ttg", "ctg", "atg");
        Genome testGto = new Genome(new File("src/test", "test.gto"));
        PegProposal prop1 = PegProposal.create(testGto, Location.create(TEST_CONTIG, "+", 1249, 1302), "hypothetical protein", 86);
        assertThat(prop1.getFunction(), equalTo("hypothetical protein"));
        assertThat(prop1.getStrength(), closeTo(0.4155, 0.0001));
        Location loc = prop1.getLoc();
        assertThat(loc.getContigId(), equalTo(TEST_CONTIG));
        assertThat(loc.getDir(), equalTo('+'));
        // Verify we found the correct endpoint.
        assertThat(loc.getEnd(), equalTo(1422));
        assertThat(loc.getRight(), equalTo(1422));
        // Verify that we contain the region.
        assertThat(loc.getLeft(), lessThanOrEqualTo(1294));
        // Verify that we extended to a start.
        String dna = testGto.getContig(TEST_CONTIG).getSequence();
        assertTrue(starts.contains(dna, loc.getBegin()));
        // Test comparison and equality.
        PegProposal prop2 = PegProposal.create(testGto, Location.create(TEST_CONTIG, "+", 1261, 1320), "serious protein", 86);
        assertThat(prop2.getFunction(), equalTo("serious protein"));
        assertThat(prop2.getStrength(), closeTo(0.5029, 0.0001));
        loc = prop2.getLoc();
        assertThat(loc.getContigId(), equalTo(TEST_CONTIG));
        assertThat(loc.getDir(), equalTo('+'));
        // Verify we found the correct endpoints.
        assertThat(loc.getEnd(), equalTo(1422));
        assertThat(loc.getRight(), equalTo(1422));
        assertThat(loc.getBegin(), equalTo(1252));
        // Test the comparisons.
        assertThat(prop2.compareTo(prop1), equalTo(0));
        assertThat(prop2, equalTo(prop1));
        assertThat(prop2.hashCode(), equalTo(prop1.hashCode()));
        assertTrue(prop1.betterThan(prop2));
        assertFalse(prop2.betterThan(prop1));
        // Test merging.
        prop1.merge(prop2);
        loc = prop1.getLoc();
        assertThat(loc.getEnd(), equalTo(1422));
        assertThat(loc.getRight(), equalTo(1422));
        assertThat(loc.getBegin(), equalTo(1252));
        assertThat(prop2, equalTo(prop1));
        assertThat(prop1.getFunction(), equalTo("serious protein"));
        assertThat(prop1.getStrength(), closeTo(0.5029, 0.0001));
        // Test a bad proposal.
        prop1 = PegProposal.create(testGto, Location.create(TEST_CONTIG, "+", 1261, 1463), "invalid protein", 0);
        assertNull(prop1);
    }

    /**
     * test proposal lists
     * @throws IOException
     */
    public void testProposalLists() throws IOException {
        Genome testGto = new Genome(new File("src/test", "test.gto"));
        PegProposalList proposals = new PegProposalList(testGto, 0.50, 80);
        // First test-- too weak
        proposals.propose(Location.create(TEST_CONTIG, "+", 1249, 1302), "long function", 69);
        assertThat(proposals.getWeakCount(), equalTo(1));
        assertThat(proposals.getProposalCount(), equalTo(0));
        assertThat(proposals.getMadeCount(), equalTo(1));
        // More evidence, will be stored
        proposals.propose(Location.create(TEST_CONTIG, "+", 1249, 1302), "long function", 138);
        assertThat(proposals.getProposalCount(), equalTo(1));
        assertThat(proposals.getMadeCount(), equalTo(2));
        // Shorter with same strength, will not be stored.
        proposals.propose(Location.create(TEST_CONTIG, "+", 1261, 1320), "short function", 114);
        assertThat(proposals.getProposalCount(), equalTo(1));
        assertThat(proposals.getMergeCount(), equalTo(0));
        assertThat(proposals.getMadeCount(), equalTo(3));
        // Shorter with more strength, gets merged.
        proposals.propose(Location.create(TEST_CONTIG, "+", 1261, 1320), "short function", 141);
        assertThat(proposals.getProposalCount(), equalTo(1));
        assertThat(proposals.getMergeCount(), equalTo(1));
        assertThat(proposals.getMadeCount(), equalTo(4));
        assertThat(proposals.getSmallCount(), equalTo(0));
        // Test min-evidence filter.
        proposals.propose(Location.create(TEST_CONTIG, "-", 100, 110), "small function", 75);
        assertThat(proposals.getProposalCount(), equalTo(1));
        assertThat(proposals.getMergeCount(), equalTo(1));
        assertThat(proposals.getMadeCount(), equalTo(5));
        assertThat(proposals.getSmallCount(), equalTo(1));
        proposals.propose(Location.create(TEST_CONTIG, "-", 100, 110), "small function", 85);
        assertThat(proposals.getProposalCount(), equalTo(2));
        assertThat(proposals.getMergeCount(), equalTo(1));
        assertThat(proposals.getMadeCount(), equalTo(6));
        assertThat(proposals.getSmallCount(), equalTo(1));
        // Add some more for variety.
        proposals.propose(Location.create(TEST_CONTIG, "+", 825851, 825853), "far protein", 163);
        proposals.propose(Location.create(TEST_CONTIG, "-", 777932, 779122), "minus protein", 600);
        proposals.propose(Location.create(TEST_CONTIG, "-", 914899, 916002), "minus 1104", 800);
        // One to reject.
        proposals.propose(Location.create(TEST_CONTIG, "+", 983222, 983349), "invalid function", 60);
        // One last weak one.
        proposals.propose(Location.create(TEST_CONTIG, "+", 905257, 905415), "weak function", 61);
        // Check the final counts.
        assertThat(proposals.getProposalCount(), equalTo(5));
        assertThat(proposals.getMergeCount(), equalTo(1));
        assertThat(proposals.getWeakCount(), equalTo(2));
        assertThat(proposals.getRejectedCount(), equalTo(1));
        assertThat(proposals.getSmallCount(), equalTo(1));
        assertThat(proposals.getMadeCount(), equalTo(11));
        // Loop through the list.
        Iterator<PegProposal> pegIter = proposals.iterator();
        assertThat(pegIter.next().getFunction(), equalTo("small function"));
        assertThat(pegIter.next().getFunction(), equalTo("short function"));
        assertThat(pegIter.next().getFunction(), equalTo("minus protein"));
        assertThat(pegIter.next().getFunction(), equalTo("far protein"));
        assertThat(pegIter.next().getFunction(), equalTo("minus 1104"));
        assertFalse(pegIter.hasNext());
    }

    /**
     * test framed location lists
     */
    public void testFramedLocations() {
        FramedLocationLists framer = new FramedLocationLists();
        // Get some locations.
        Location loc1p = Location.create("c1", "+", 100, 200);
        Location loc2p = Location.create("c1", "+", 103, 203);
        Location loc3p = Location.create("c1", "+", 101, 201);
        Location loc4p = Location.create("c1", "+", 104, 204);
        Location loc5p = Location.create("c1", "+", 102, 202);
        Location loc6p = Location.create("c1", "+", 105, 205);
        Location loc1m = Location.create("c1", "-", 100, 200);
        Location loc2m = Location.create("c1", "-", 103, 203);
        Location loc3m = Location.create("c1", "-", 101, 201);
        Location loc4m = Location.create("c1", "-", 104, 204);
        Location loc5m = Location.create("c1", "-", 102, 202);
        Location loc6m = Location.create("c1", "-", 105, 205);
        // These sets are used to determine which locations go with with target.
        HashSet<Location> target1 = new HashSet<Location>();
        target1.add(loc1p); target1.add(loc2p); target1.add(loc4p); target1.add(loc1m); target1.add(loc2m); target1.add(loc5m);
        HashSet<Location> target2 = new HashSet<Location>();
        target2.add(loc3p); target2.add(loc5p); target2.add(loc6p); target2.add(loc3m); target2.add(loc4m); target2.add(loc6m);
        // Now put these in the framer.
        for (Location loc : target1) framer.connect("t1", loc);
        for (Location loc : target2) framer.connect("t2", loc);
        assertThat(framer.size(), equalTo(12));
        // This set is used to track what we found.
        Collection<Location> found = new ArrayList<Location>();
        // Iterate through the framer.
        for (FramedLocationLists.Report report : framer) {
            SortedLocationList locs = report.getList();
            // Insure we got a valid target back.
            assertThat(report.getId(), anyOf(equalTo("t1"), equalTo("t2")));
            // Insure all the locations are for this target and belong to the same frame.
            Set<Location> target = (report.getId().contentEquals("t1") ? target1 : target2);
            Frame frm = Frame.XX;
            for (Location loc : locs) {
                found.add(loc);
                assertThat(target, hasItem(loc));
                if (frm == Frame.XX)
                    frm = loc.getFrame();
                else
                    assertThat(loc.getFrame(), equalTo(frm));
            }
        }
        // Insure we found all the entries and only those entries.
        assertThat(found.size(), equalTo(12));
        for (Location loc : target1) assertThat(found, hasItem(loc));
        for (Location loc : target2) assertThat(found, hasItem(loc));
        // Clear the lists.
        framer.clear();
        Iterator<FramedLocationLists.Report> iter = framer.iterator();
        assertFalse(iter.hasNext());
        // Insert just target 1.
        for (Location loc : target1) framer.connect("t1", loc);
        // Iterate again, insuring we only have target 1.
        for (FramedLocationLists.Report report : framer) {
            assertThat(report.getId(), equalTo("t1"));
            for (Location loc : report.getList()) {
                assertThat(target1, hasItem(loc));
            }
        }
        assertThat(framer.size(), equalTo(6));
    }
}
