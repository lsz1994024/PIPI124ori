package ProteomicsLibrary;

import com.google.common.collect.Multimap;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DbTool {

    private static final Random random = new Random();

    private Map<String, String> proteinSequenceMap = new HashMap<>();
    private Map<String, String> proteinAnnotateMap = new HashMap<>();

    public DbTool(String dbName, String databaseType) throws IOException {
        String id = "";
        String annotate;
        StringBuilder sequence = new StringBuilder(99999);
        databaseType = databaseType.trim().toLowerCase();

        boolean newPro = true;

        Pattern headerPattern;
        if (databaseType.contentEquals("tair")) {
            headerPattern = Pattern.compile("^>([^\\s]+)[\\s|]*(.*)$");
        } else if (databaseType.contentEquals("uniprot") || databaseType.contentEquals("swissprot")) {
            headerPattern = Pattern.compile("^>([^ ]+) *(.*)$");
        } else if (databaseType.contentEquals("nextprot")) {
            headerPattern = Pattern.compile("^>([^ ]+) *(.*)");
        } else if (databaseType.contentEquals("contaminants") || databaseType.contentEquals("itag") || databaseType.contentEquals("refseq")) {
            headerPattern = Pattern.compile("^>([^ ]+) *(.*)$");
        } else if (databaseType.contentEquals("others")) {
            headerPattern = Pattern.compile("^>(.+)$");
        } else {
            throw new NullPointerException(String.format(Locale.US, "Incorrect database type (%s) in the parameter file.", databaseType));
        }

        BufferedReader dbReader;
        if (databaseType.contentEquals("contaminants")) {
            InputStream inputStream = getClass().getClassLoader().getResourceAsStream("contaminants.fasta");
            dbReader = new BufferedReader(new InputStreamReader(inputStream));
        } else {
            dbReader = new BufferedReader(new FileReader(dbName));
        }
        String line;
        while ((line = dbReader.readLine()) != null) {
            line = line.trim();
            Matcher headMatcher = headerPattern.matcher(line);
            if (headMatcher.matches()) {

                if (!newPro) {
                    proteinSequenceMap.put(id, sequence.toString());
                }
                id = headMatcher.group(1).trim();
                if (databaseType.contentEquals("others")) {
                    annotate = id;
                } else {
                    annotate = headMatcher.group(2).trim();
                }
                proteinAnnotateMap.put(id, annotate);
                newPro = true;
            } else if (!line.isEmpty()) {

                if (newPro) {
                    sequence = new StringBuilder(99999);
                    sequence.append(line);
                    newPro = false;
                } else {
                    sequence.append(line);
                }
            }
        }
        dbReader.close();

        proteinSequenceMap.put(id, sequence.toString());
    }

    public Map<String, String> getProteinSequenceMap() {
        return proteinSequenceMap;
    }

    public Map<String, String> getProteinAnnotateMap() {
        return proteinAnnotateMap;
    }

    public static String shuffleSeq(String sequence, String cleavageSite, String protectionSite, boolean cleavageFromCTerm) {


        String sequenceToBeShuffled;
        if (sequence.startsWith("M")) {
            sequenceToBeShuffled = sequence.substring(1);
        } else {
            sequenceToBeShuffled = sequence;
        }

        Pattern digestSitePattern = MassTool.getDigestSitePattern(cleavageSite, protectionSite, cleavageFromCTerm);
        Set<Integer> cutSiteSet = new HashSet<>();
        Matcher matcher = digestSitePattern.matcher(sequenceToBeShuffled);
        while (matcher.find()) {
            cutSiteSet.add(matcher.start());
        }
        char[] tempArray = sequenceToBeShuffled.toCharArray();
        int idx = 0;
        while (idx < tempArray.length - 1) {
            if (!cutSiteSet.contains(idx) && !cutSiteSet.contains(idx + 1)) {
                char temp = tempArray[idx];
                tempArray[idx] = tempArray[idx + 1];
                tempArray[idx + 1] = temp;
                idx += 2;
            } else {
                ++idx;
            }
        }

        if (sequence.startsWith("M")) {
            return "M" + String.valueOf(tempArray);
        } else {
            return String.valueOf(tempArray);
        }
    }

    public static Character[] getLeftRightFlank(String peptide, Multimap<String, String> peptideProteinMap, Map<String, String> proteinSequenceMap, String cleavageSite, String protectionSite, boolean cleavageFromCTerm) throws Exception {
        Character[] leftRightFlank = new Character[2];
        String peptideString = DbTool.getSequenceOnly(peptide);
        for (String proteinId : peptideProteinMap.get(peptide)) {
            String proteinSequence = proteinSequenceMap.get(proteinId);
            int startIdx = proteinSequence.indexOf(peptideString);
            while (startIdx >= 0) {
                if (startIdx == 0 || ((startIdx == 1 && proteinSequence.charAt(0) == 'M'))) {
                    int tempIdx = startIdx + peptideString.length();
                    if (tempIdx < proteinSequence.length()) {
                        leftRightFlank[1] = proteinSequence.charAt(tempIdx);
                        if ((cleavageFromCTerm && !protectionSite.contains(leftRightFlank[1].toString())) || (!cleavageFromCTerm && cleavageSite.contains(leftRightFlank[1].toString()))) {
                            leftRightFlank[0] = '-';
                            break;
                        } else {
                            leftRightFlank[1] = null;
                        }
                    } else if (tempIdx == proteinSequence.length()) {
                        leftRightFlank[0] = '-';
                        leftRightFlank[1] = '-';
                        break;
                    } else {
                        throw new Exception(String.format(Locale.US, "The peptide %s is longer than its protein %s.", peptideString, proteinSequence));
                    }
                } else if (startIdx == proteinSequence.length() - peptideString.length()) {
                    leftRightFlank[0] = proteinSequence.charAt(startIdx - 1);
                    if ((cleavageFromCTerm && cleavageSite.contains(leftRightFlank[0].toString())) || (!cleavageFromCTerm && !protectionSite.contains(leftRightFlank[0].toString()))) {
                        leftRightFlank[1] = '-';
                        break;
                    } else {
                        leftRightFlank[0] = null;
                    }
                } else {
                    leftRightFlank[0] = proteinSequence.charAt(startIdx - 1);
                    leftRightFlank[1] = proteinSequence.charAt(startIdx + peptideString.length());
                    if ((cleavageFromCTerm && cleavageSite.contains(leftRightFlank[0].toString()) && !protectionSite.contains(leftRightFlank[1].toString())) || (!cleavageFromCTerm && cleavageSite.contains(leftRightFlank[1].toString()) && !protectionSite.contains(leftRightFlank[0].toString()))) {
                        break;
                    } else {
                        leftRightFlank[0] = null;
                        leftRightFlank[1] = null;
                    }
                }
                startIdx = proteinSequence.indexOf(peptideString, startIdx + 1);
            }

            if (leftRightFlank[0] != null && leftRightFlank[1] != null) {
                return leftRightFlank;
            }
        }
        return null;
    }

    public static Set<Integer> findPeptideLocation(String proteinSequence, String peptide, String cutSite, String protectSite) throws NullPointerException {
        peptide = getSequenceOnly(peptide.trim());
        Set<Integer> output = new HashSet<>();
        int idx = proteinSequence.indexOf(peptide);
        while (idx >= 0) {
            if ((idx == 0 || cutSite.contains(proteinSequence.substring(idx - 1, idx)) || (idx == 1 && proteinSequence.charAt(0) == 'M')) && (idx + peptide.length() == proteinSequence.length() || !protectSite.contains(proteinSequence.substring(idx + peptide.length(), idx + peptide.length() + 1)))) {
                output.add(idx);
            }
            idx = proteinSequence.indexOf(peptide, idx + 1);
        }
        return output;
    }

    public static String getPtmFreePeptide(String peptide) {
        return peptide.replaceAll("[^A-Znc]+", "");
    }

    public static String getSequenceOnly(String peptide) {
        return peptide.replaceAll("[^A-Z]+", "");
    }
}
