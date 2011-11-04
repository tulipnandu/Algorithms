/**
 *
 */
package assign;

import assign.Consensus;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;


/**
 * @author Yogesh
 *
 */
@SuppressWarnings("deprecation")
public class GreedyMotifSearch2 {
    static SequenceDB seqDB;
    static int length;
    static String inFile;
    static List<Motif> motifs;
    private static List<String> allSeq;
    static List<Consensus> cons;
    private static final String[] aa = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" };

    /**
     * @param args
     * @throws BioException
     * @throws FileNotFoundException
     */
    public static void main(String[] args) throws FileNotFoundException, BioException {
        // TODO Auto-generated method stub
        long time = System.currentTimeMillis();
        GreedyMotifSearch2 gms2 = new GreedyMotifSearch2();

        //inFile = "./resources/better_sample.fasta";
        //inFile = "./resources/zinc_finger.fasta";
        //inFile = "./resources/short.fasta";
        if (args.length == 2) {
            inFile = args [0];
            length = Integer.parseInt(args [1]);

            readInput("PROTEIN");

            Motif motif = gms2.searchMotif();

            //System.out.println(motif.getMotifString() + "\t" + motif.getScore());
            List<Consensus> newCons = sortConsensusMotifs();
            String res = print(newCons, 20);
            System.out.println(res);

            time = System.currentTimeMillis() - time;
            System.out.println();
            System.out.println("Approx. execution time for GreedyMotifSearch: " + (time / 1000) + " seconds (" + time + " milliseconds)");
        } else {
            System.out.println("Usage: \njava -jar GreedyMotifSearch.jar <input.fasta> <l-mer>");
        }

        //length = 25;
    }

    /**
     *
     * @param type
     * @throws BioException
     * @throws FileNotFoundException
     */
    private static void readInput(String type) throws BioException, FileNotFoundException {
        BufferedInputStream bis = new BufferedInputStream(new FileInputStream(inFile));

        Alphabet alpha = AlphabetManager.alphabetForName(type);

        seqDB = SeqIOTools.readFasta(bis, alpha);
    }

    /**
     *
     * @throws NoSuchElementException
     * @throws BioException
     */
    private static void createFunction() throws NoSuchElementException, BioException {
        motifs = new ArrayList<Motif>();

        int bmScore = 0;
        int bStart1 = 0;
        int bStart2 = 0;
        String bM1 = "";
        String bM2 = "";

        allSeq = new ArrayList<String>();

        SequenceIterator seqIter = seqDB.sequenceIterator();

        while (seqIter.hasNext()) {
            Sequence sq = seqIter.nextSequence();

            String seq = sq.seqString();
            //String name = sq.getName();
            //System.out.println(name + "\t" + seq);
            allSeq.add(seq);
        }

        String seq1 = allSeq.get(0);

        String seq2 = allSeq.get(1);

        //System.out.println(seq1 + "\n" + seq2);
        if (length < seq1.length()) {
            for (int i = 0; i < (seq1.length() - length + 1); i++) {
                String s1 = seq1.substring(i, i + length);
                int start1 = i;

                if (length < seq2.length()) {
                    for (int j = 0; j < (seq2.length() - length + 1); j++) {
                        String s2 = seq2.substring(j, j + length);
                        int start2 = j;
                        int dH = getHammingDistance(s1, s2);

                        //System.out.println(s1 + "\t" + s2 + "\t" + dH);
                        if ((length - dH) > bmScore) {
                            bmScore = (length - dH);
                            bM1 = s1;
                            bM2 = s2;
                            bStart1 = start1 + 1;
                            bStart2 = start2 + 1;
                        }
                    }
                } else {
                    System.out.println("Motif cannot be longer than the Sequence!!!");

                    break;
                }
            }
        } else {
            System.out.println("Motif cannot be longer than the Sequence!!!");
        }

        //System.out.println("Maximization function created");
        //System.out.println("Following 2 are the seeds");
        motifs.add(new Motif(bM1, bStart1, bStart1 + length, bmScore));
        motifs.add(new Motif(bM2, bStart2, bStart2 + length, bmScore));

        //System.out.println(bM1 + ", " + bStart1 + "-" + (bStart1 + length) + ", " + bmScore);
        //System.out.println(bM2 + ", " + bStart2 + "-" + (bStart2 + length) + ", " + bmScore);
        //System.out.println();

        //System.out.println (motifs.size());
    }

    private static List<Consensus> sortConsensusMotifs() {
        /*List<Consensus> newCons = new ArrayList<Consensus>();
        
        for (int c1 = 1; c1 < cons.size(); c1++) {
            for (int c2 = 0; c2 < (cons.size() - 1); c2++) {
                Consensus tC1 = cons.get(c1);
                Consensus tC2 = cons.get(c2);
            }
        }
        
        return newCons;*/

        //System.out.println("Sorting the Consensus List");
        Collections.sort(cons,
            new Comparator<Consensus>() {
                public int compare(Consensus c1, Consensus c2) {
                    int a = c1.getConsensusScore();
                    int b = c2.getConsensusScore();

                    if (a < b) {
                        return 1;
                    } else if (a > b) {
                        return -1;
                    } else {
                        return 0;
                    }
                }
            });

        return cons;
    }

    /**
     *
     * @return Motif <code>String</code> obtained from Greedy search
     * @throws NoSuchElementException
     * @throws BioException
     */
    public Motif searchMotif() throws NoSuchElementException, BioException {
        createFunction();

        int SCORE = 0;
        String MOTIF = "";
        int START = 0;
        int END = 0;

        cons = new ArrayList<Consensus>();

        for (int i = 2; i < allSeq.size(); i++) {
            String seq = allSeq.get(i);

            //System.out.println("Sequence #" + (i + 1));
            if (seq.length() > length) {
                for (int j = 0; j < (seq.length() - length + 1); j++) {
                    String s = seq.substring(j, j + length);

                    String[] consensus = getScore(s);
                    int bestScore = Integer.parseInt(consensus [1]);

                    /*String[] str1 = new String[motifs.size()+1];
                    int n = 0;
                    for (int m = 0; m < motifs.size(); m++) {
                        str1 [m] = motifs.get(m).getMotifString();
                        n = m;
                    }
                    
                    str1 [n+1] = s;
                    
                    String[] best = getBestScore(str1, (motifs.size()+1), length);
                    int bestScore = Integer.parseInt(best [1]);*/
                    if (i == (allSeq.size() - 1)) {
                        motifs.add(new Motif(consensus [0], j + 1, (j + 1 + length), bestScore));

                        Motif m = getConsensusMotif();
                        String mS = m.getMotifString();
                        int mSco = m.getScore();
                        int[] start = new int[motifs.size() + 1];

                        for (int st = 0; st < motifs.size(); st++) {
                            int strt = motifs.get(st).getStartPosition();
                            start [st] = strt;
                        }

                        cons.add(new Consensus(mS, start, mSco));
                        motifs.remove(motifs.size() - 1);
                    }

                    //System.out.println(s + "-" + score);
                    if (bestScore > SCORE) {
                        SCORE = bestScore;
                        MOTIF = s;
                        START = j + 1;
                        END = START + length;
                    }
                }
            } else {
                //System.out.println("Motif cannot be longer than the Sequence!!!\nNot considering this sequence");
                continue;
            }

            motifs.add(new Motif(MOTIF, START, END, SCORE));

            //System.out.println(MOTIF + ", " + START + "-" + END + ", " + SCORE);

            //MOTIF = "";
            //SCORE = 0;
            //START = 0;
            //END = 0;
        }

        System.out.println();

        //System.out.println("No. of Motifs in Profile:  " + motifs.size());
        return getConsensusMotif();
    }

    /**
     *
     * @param s
     * @return <code>int</code> for score from Profile
     */
    public String[] getScore(String s) {
        List<Motif> tempMotifs = new ArrayList<Motif>();

        for (int m = 0; m < motifs.size(); m++) {
            tempMotifs.add(motifs.get(m));
        }

        //tempMotifs = motifs;
        tempMotifs.add(new Motif(s, 0, 0, 0));

        //System.out.println("No. of Sequences in Orig: " + motifs.size());
        //System.out.println ("No. of Sequences in Temp: "+tempMotifs.size());
        String[] consensus = getAlphabetCount(tempMotifs);

        //int[] aa = new int[20];
        return consensus;
    }

    /**
     *
     * @param aa
     * @param aaFreq
     * @return
     */
    private static String[] getConsensusAA(Map<String, Integer> aaFreq) {
        String[] bestAA = new String[2];
        int high = 0;

        for (int a = 0; a < aa.length; a++) {
            String amino = aa [a];
            int val = aaFreq.get(aa [a]);

            //System.out.println (amino+"\t"+val);
            if (val > high) {
                bestAA [0] = amino;
                bestAA [1] = String.valueOf(val);
            }
        }

        //System.out.println(bestAA [0] + " (" + bestAA [1] + ")");

        //System.out.print(",");
        return bestAA;
    }

    /**
     *
     * @return Consensus motif as <code>String</code>
     */
    public Motif getConsensusMotif() {
        Motif motif;
        String[] consensus = getAlphabetCount(motifs);
        motif = new Motif(consensus [0], 0, 0, Integer.parseInt(consensus [1]));

        return motif;
    }

    /**
     *
     * @param s1
     * @param s2
     * @return Integer value of Hamming distance
     */
    private static int getHammingDistance(String s1, String s2) {
        if (s1.length() != s2.length()) {
            return -1;
        }

        int counter = 0;

        for (int h = 0; h < s1.length(); h++) {
            if (s1.charAt(h) != s2.charAt(h)) {
                counter++;
            }
        }

        return counter;
    }

    public static String print(List<Consensus> nCons, int n) {
        StringBuffer sb = new StringBuffer();

        for (int p = 0; p < n; p++) {
            String toPrint = "";
            Consensus con = nCons.get(p);
            int[] start = con.getStartPositions();
            toPrint = "[";

            for (int s = 0; s < start.length; s++) {
                toPrint = toPrint + " " + start [s];
            }

            toPrint = toPrint + " ]";
            toPrint = toPrint + "\t" + con.getConsensusScore();
            toPrint = toPrint + "\t" + con.getMotifString();
            //System.out.println (toPrint);
            sb.append(toPrint + "\n");

            //System.out.println (sb.toString());
        }

        return sb.toString();
    }

    /**
     *
     * @param tMotifs
     * @param getAA
     * @return
     */
    private static String[] getAlphabetCount(List<Motif> tMotifs) {
        Map<String, Integer> col = new HashMap<String, Integer>();
        int score = 0;
        String consensus = "";
        String[] mot = new String[2];

        for (int c1 = 0; c1 < length; c1++) {
            for (int a = 0; a < aa.length; a++) {
                col.put(aa [a], 0);
            }

            for (int c2 = 0; c2 < tMotifs.size(); c2++) {
                String tM = tMotifs.get(c2).getMotifString();
                String ch = String.valueOf(tM.charAt(c1));

                if (ch.equals("A")) {
                    String key = "A";
                    int A = col.get(key);
                    A = A + 1;
                    col.remove("A");
                    col.put("A", A);
                } else if (ch.equals("C")) {
                    String key = "C";
                    int C = col.get(key);
                    C = C + 1;
                    col.remove("C");
                    col.put("C", C);
                } else if (ch.equals("D")) {
                    String key = "D";
                    int D = col.get(key);
                    D = D + 1;
                    col.remove("D");
                    col.put("D", D);
                } else if (ch.equals("E")) {
                    String key = "E";
                    int E = col.get(key);
                    E = E + 1;
                    col.remove("E");
                    col.put("E", E);
                } else if (ch.equals("F")) {
                    String key = "F";
                    int F = col.get(key);
                    F = F + 1;
                    col.remove("F");
                    col.put("F", F);
                } else if (ch.equals("G")) {
                    String key = "G";
                    int G = col.get(key);
                    G = G + 1;
                    col.remove("G");
                    col.put("G", G);
                } else if (ch.equals("H")) {
                    String key = "H";
                    int H = col.get(key);
                    H = H + 1;
                    col.remove("H");
                    col.put("H", H);
                } else if (ch.equals("I")) {
                    String key = "I";
                    int I = col.get(key);
                    I = I + 1;
                    col.remove("I");
                    col.put("I", I);
                } else if (ch.equals("J")) {
                    String key = "J";
                    int J = col.get(key);
                    J = J + 1;
                    col.remove("J");
                    col.put("J", J);
                } else if (ch.equals("K")) {
                    String key = "K";
                    int K = col.get(key);
                    K = K + 1;
                    col.remove("K");
                    col.put("K", K);
                } else if (ch.equals("L")) {
                    String key = "L";
                    int L = col.get(key);
                    L = L + 1;
                    col.remove("L");
                    col.put("L", L);
                } else if (ch.equals("M")) {
                    String key = "M";
                    int M = col.get(key);
                    M = M + 1;
                    col.remove("M");
                    col.put("M", M);
                } else if (ch.equals("N")) {
                    String key = "N";
                    int N = col.get(key);
                    N = N + 1;
                    col.remove("N");
                    col.put("N", N);
                } else if (ch.equals("P")) {
                    String key = "P";
                    int P = col.get(key);
                    P = P + 1;
                    col.remove("P");
                    col.put("P", P);
                } else if (ch.equals("Q")) {
                    String key = "Q";
                    int Q = col.get(key);
                    Q = Q + 1;
                    col.remove("Q");
                    col.put("Q", Q);
                } else if (ch.equals("R")) {
                    String key = "R";
                    int R = col.get(key);
                    R = R + 1;
                    col.remove("R");
                    col.put("R", R);
                } else if (ch.equals("S")) {
                    String key = "S";
                    int S = col.get(key);
                    S = S + 1;
                    col.remove("S");
                    col.put("S", S);
                } else if (ch.equals("T")) {
                    String key = "T";
                    int T = col.get(key);
                    T = T + 1;
                    col.remove("T");
                    col.put("T", T);
                } else if (ch.equals("V")) {
                    String key = "V";
                    int V = col.get(key);
                    V = V + 1;
                    col.remove("V");
                    col.put("V", V);
                } else if (ch.equals("W")) {
                    String key = "W";
                    int W = col.get(key);
                    W = W + 1;
                    col.remove("W");
                    col.put("W", W);
                } else if (ch.equals("Y")) {
                    String key = "Y";
                    int Y = col.get(key);
                    Y = Y + 1;
                    col.remove("Y");
                    col.put("Y", Y);
                }
            }

            String[] bestAA = getConsensusAA(col);
            consensus = consensus + bestAA [0];
            //System.out.println (bestAA[1]);
            score = score + Integer.valueOf(bestAA [1]);
        }

        mot [0] = consensus;
        mot [1] = Integer.toString(score);

        return mot;
    }

    /*int A = 0;
    int C = 0;
    int D = 0;
    int E = 0;
    int F = 0;
    int G = 0;
    int H = 0;
    int I = 0;
    int J = 0;
    int K = 0;
    int L = 0;
    int M = 0;
    int N = 0;
    int P = 0;
    int Q = 0;
    int R = 0;
    int S = 0;
    int T = 0;
    int V = 0;
    int W = 0;
    int Y = 0;*/
    public static String[] getBestScore(String[] str1, int row, int col) {
        String s1;

        //    DataInputStream dis = new DataInputStream(System.in);
        int i = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int j = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int //            "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
         Acount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Ccount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Dcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Ecount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Fcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Gcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Hcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Icount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Kcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Lcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Mcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Ncount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Pcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Qcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Rcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Scount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Tcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Vcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Wcount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int Ycount = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int max = 0;

        //    DataInputStream dis = new DataInputStream(System.in);
        int finalscore = 0;
        int[] A = new int[col + 1];
        int[] C = new int[col + 1];
        int[] D = new int[col + 1];
        int[] E = new int[col + 1];
        int[] F = new int[col + 1];
        int[] G = new int[col + 1];
        int[] H = new int[col + 1];
        int[] I = new int[col + 1];
        int[] K = new int[col + 1];
        int[] L = new int[col + 1];
        int[] M = new int[col + 1];
        int[] N = new int[col + 1];
        int[] P = new int[col + 1];
        int[] Q = new int[col + 1];
        int[] R = new int[col + 1];
        int[] S = new int[col + 1];
        int[] T = new int[col + 1];
        int[] V = new int[col + 1];
        int[] W = new int[col + 1];
        int[] Y = new int[col + 1];

        int[] scorearr = new int[col + 1];
        char[] str = new char[500];
        char[] scoreseq = new char[col + 1];
        char[][] arr = new char[500][500];

        for (i = 0; i < row; i++) {
            s1 = str1 [i];
            str = s1.toCharArray();

            for (j = 0; j < col; j++) {
                arr [i] [j] = str [j];
            }
        }

        for (i = 0; i < row; i++) {
            for (j = 0; j < col; j++) {
                //System.out.print(arr [i] [j] + "\t");
            }

            //System.out.println("\n");
        }

        for (i = 0; i <= row; i++) {
            for (j = 0; j < (col + 1); j++) {
                if ((arr [j] [i] == 'A') || (arr [j] [i] == 'a')) {
                    Acount = Acount + 1;
                } else if ((arr [j] [i] == 'C') || (arr [j] [i] == 'c')) {
                    Ccount = Ccount + 1;
                } else if ((arr [j] [i] == 'D') || (arr [j] [i] == 'd')) {
                    Dcount = Dcount + 1;
                } else if ((arr [j] [i] == 'E') || (arr [j] [i] == 'e')) {
                    Ecount = Ecount + 1;
                } else if ((arr [j] [i] == 'F') || (arr [j] [i] == 'f')) {
                    Fcount = Fcount + 1;
                } else if ((arr [j] [i] == 'G') || (arr [j] [i] == 'g')) {
                    Gcount = Gcount + 1;
                } else if ((arr [j] [i] == 'H') || (arr [j] [i] == 'h')) {
                    Hcount = Hcount + 1;
                } else if ((arr [j] [i] == 'I') || (arr [j] [i] == 'i')) {
                    Icount = Icount + 1;
                } else if ((arr [j] [i] == 'K') || (arr [j] [i] == 'k')) {
                    Kcount = Kcount + 1;
                } else if ((arr [j] [i] == 'L') || (arr [j] [i] == 'l')) {
                    Lcount = Lcount + 1;
                } else if ((arr [j] [i] == 'M') || (arr [j] [i] == 'm')) {
                    Mcount = Mcount + 1;
                } else if ((arr [j] [i] == 'N') || (arr [j] [i] == 'n')) {
                    Ncount = Ncount + 1;
                } else if ((arr [j] [i] == 'P') || (arr [j] [i] == 'p')) {
                    Pcount = Pcount + 1;
                } else if ((arr [j] [i] == 'Q') || (arr [j] [i] == 'q')) {
                    Qcount = Qcount + 1;
                } else if ((arr [j] [i] == 'R') || (arr [j] [i] == 'r')) {
                    Rcount = Rcount + 1;
                } else if ((arr [j] [i] == 'S') || (arr [j] [i] == 's')) {
                    Scount = Scount + 1;
                } else if ((arr [j] [i] == 'T') || (arr [j] [i] == 't')) {
                    Tcount = Tcount + 1;
                } else if ((arr [j] [i] == 'V') || (arr [j] [i] == 'v')) {
                    Vcount = Vcount + 1;
                } else if ((arr [j] [i] == 'W') || (arr [j] [i] == 'w')) {
                    Wcount = Wcount + 1;
                } else if ((arr [j] [i] == 'Y') || (arr [j] [i] == 'y')) {
                    Ycount = Ycount + 1;
                }
            }

            A [i + 1] = Acount;
            Acount = 0;
            C [i + 1] = Ccount;
            Ccount = 0;
            D [i + 1] = Dcount;
            Dcount = 0;
            E [i + 1] = Ecount;
            Ecount = 0;
            F [i + 1] = Fcount;
            Fcount = 0;
            G [i + 1] = Gcount;
            Gcount = 0;
            H [i + 1] = Hcount;
            Hcount = 0;
            I [i + 1] = Icount;
            Icount = 0;
            K [i + 1] = Kcount;
            Kcount = 0;
            L [i + 1] = Lcount;
            Lcount = 0;
            M [i + 1] = Mcount;
            Mcount = 0;
            N [i + 1] = Ncount;
            Ncount = 0;
            P [i + 1] = Pcount;
            Pcount = 0;
            Q [i + 1] = Qcount;
            Qcount = 0;
            R [i + 1] = Rcount;
            Rcount = 0;
            S [i + 1] = Scount;
            Scount = 0;
            T [i + 1] = Tcount;
            Tcount = 0;
            V [i + 1] = Vcount;
            Vcount = 0;
            W [i + 1] = Wcount;
            Wcount = 0;
            Y [i + 1] = Ycount;
            Ycount = 0;
        }

        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(65 + "\t");
            } else {
                //System.out.print(A [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(67 + "\t");
            } else {
                //System.out.print(C [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(68 + "\t");
            } else {
                //System.out.print(D [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(69 + "\t");
            } else {
                //System.out.print(E [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(70 + "\t");
            } else {
                //System.out.print(F [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(71 + "\t");
            } else {
                //System.out.print(G [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(72 + "\t");
            } else {
                //System.out.print(H [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(73 + "\t");
            } else {
                //System.out.print(I [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(75 + "\t");
            } else {
                //System.out.print(K [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(76 + "\t");
            } else {
                //System.out.print(L [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(77 + "\t");
            } else {
                //System.out.print(M [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(78 + "\t");
            } else {
                //System.out.print(N [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(80 + "\t");
            } else {
                //System.out.print(P [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(81 + "\t");
            } else {
                //System.out.print(Q [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(82 + "\t");
            } else {
                //System.out.print(R [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(83 + "\t");
            } else {
                //System.out.print(S [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(84 + "\t");
            } else {
                // System.out.print(T [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(86 + "\t");
            } else {
                //System.out.print(V [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(87 + "\t");
            } else {
                //System.out.print(W [i] + "\t");
            }
        }

        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            if (i == 0) {
                //System.out.print(89 + "\t");
            } else {
                //System.out.print(Y [i] + "\t");
            }
        }

        //System.out.print("\n");
        //System.out.print("\n");
        for (i = 0; i < (col + 1); i++) {
            max = A [i];
            scoreseq [i] = 'A';

            if (C [i] > max) {
                max = C [i];
                scoreseq [i] = 'C';
            }

            if (D [i] > max) {
                max = D [i];
                scoreseq [i] = 'D';
            }

            if (E [i] > max) {
                max = E [i];
                scoreseq [i] = 'E';
            }

            if (F [i] > max) {
                max = F [i];
                scoreseq [i] = 'F';
            }

            if (G [i] > max) {
                max = G [i];
                scoreseq [i] = 'G';
            }

            if (H [i] > max) {
                max = H [i];
                scoreseq [i] = 'H';
            }

            if (I [i] > max) {
                max = I [i];
                scoreseq [i] = 'I';
            }

            if (K [i] > max) {
                max = K [i];
                scoreseq [i] = 'K';
            }

            if (L [i] > max) {
                max = L [i];
                scoreseq [i] = 'L';
            }

            if (M [i] > max) {
                max = M [i];
                scoreseq [i] = 'M';
            }

            if (N [i] > max) {
                max = N [i];
                scoreseq [i] = 'N';
            }

            if (P [i] > max) {
                max = P [i];
                scoreseq [i] = 'P';
            }

            if (Q [i] > max) {
                max = Q [i];
                scoreseq [i] = 'Q';
            }

            if (R [i] > max) {
                max = R [i];
                scoreseq [i] = 'R';
            }

            if (S [i] > max) {
                max = S [i];
                scoreseq [i] = 'S';
            }

            if (T [i] > max) {
                max = T [i];
                scoreseq [i] = 'T';
            }

            if (V [i] > max) {
                max = V [i];
                scoreseq [i] = 'V';
            }

            if (W [i] > max) {
                max = W [i];
                scoreseq [i] = 'W';
            }

            if (Y [i] > max) {
                max = Y [i];
                scoreseq [i] = 'Y';
            }

            scorearr [i] = max;
            max = 0;
        }

        for (i = 1; i < (col + 1); i++) {
            finalscore = finalscore + scorearr [i];

            //System.out.print(scoreseq [i] + "\t");
        }

        //System.out.print("\n Score= " + finalscore);
        String consensus = new String(scoreseq);
        String[] scoring = new String[2];
        scoring [0] = consensus;
        scoring [1] = Integer.toString(finalscore);

        //System.out.print(scoring [0]);
        //System.out.print(scoring [1]);
        return (scoring);
    }
}
