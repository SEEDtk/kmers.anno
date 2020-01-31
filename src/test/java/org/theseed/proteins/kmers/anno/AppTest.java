package org.theseed.proteins.kmers.anno;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.apache.commons.lang3.StringUtils;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.locations.PegProposal;
import org.theseed.proteins.CodonSet;
import org.theseed.proteins.DnaTranslator;
import org.theseed.proteins.kmers.KmerReference;

/**
 * Unit test for simple App.
 */
public class AppTest extends TestCase
{
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
        CountMap<KmerReference> kmerCounts = KmerReference.countContigKmers(smallGto);
        // Verify that we found the kmers in all the right places.
        for (CountMap<KmerReference>.Count kmerCount : kmerCounts.counts()) {
            KmerReference kmerRef = kmerCount.getKey();
            String kmer = kmerRef.getKmer();
            Location loc = kmerRef.getLoc();
            assertFalse("Invalid characters in " + kmerRef, StringUtils.containsAny(kmer, 'X', '*'));
            assertThat("Kmer does not match location in " + kmerRef, kmer, equalTo(xlator.translate(smallGto.getDna(loc))));
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
        PegProposal prop1 = PegProposal.create(testGto, Location.create("51203.13.con.0001", "+", 1294, 1302), "hypothetical protein", 100.0);
        assertThat(prop1.getFunction(), equalTo("hypothetical protein"));
        assertThat(prop1.getStrength(), equalTo(100.0));
        Location loc = prop1.getLoc();
        assertThat(loc.getContigId(), equalTo("51203.13.con.0001"));
        assertThat(loc.getDir(), equalTo('+'));
        // Verify we found the correct endpoint.
        assertThat(loc.getEnd(), equalTo(1422));
        assertThat(loc.getRight(), equalTo(1422));
        // Verify that we contain the region.
        assertThat(loc.getLeft(), lessThanOrEqualTo(1294));
        // Verify that we extended to a start.
        String dna = testGto.getContig("51203.13.con.0001").getSequence();
        assertTrue(starts.contains(dna, loc.getBegin()));
        // Test comparison and equality.
        PegProposal prop2 = PegProposal.create(testGto, Location.create("51203.13.con.0001", "+", 1261, 1320), "serious protein", 80.0);
        assertThat(prop2.getFunction(), equalTo("serious protein"));
        assertThat(prop2.getStrength(), equalTo(80.0));
        loc = prop2.getLoc();
        assertThat(loc.getContigId(), equalTo("51203.13.con.0001"));
        assertThat(loc.getDir(), equalTo('+'));
        // Verify we found the correct endpoints.
        assertThat(loc.getEnd(), equalTo(1422));
        assertThat(loc.getRight(), equalTo(1422));
        assertThat(loc.getBegin(), equalTo(1252));
        // Test the comparisons.
        assertThat(prop2.compareTo(prop1), equalTo(0));
        assertThat(prop2, equalTo(prop1));
        assertThat(prop2.hashCode(), equalTo(prop1.hashCode()));
        // Test a bad proposal.
        prop1 = PegProposal.create(testGto, Location.create("51203.13.con.0001", "+", 1261, 1463), "invalid protein", 0.0);
        assertNull(prop1);
    }



}
