import java.io.*;
import java.util.*;


public class TextMasterEM {
	private double [][] _pZgivenW;
	private double [][] _pCgivenZ;
	public final int WORDS = 122;
	public final int CONTEXTS = 23804;
	public int TOPICS;

	public TextMasterEM(int topics) {
		TOPICS = topics;
		_pZgivenW = new double[WORDS][TOPICS];
		_pCgivenZ = new double[TOPICS][CONTEXTS];
	}
	
	//returns keys in ascending order of values
	//allows duplicate values (defaults to input order)
	public static ArrayList<String> sortKeysByVal(HashMap<String, Double> hm) {
		ArrayList<String> out = new ArrayList<String>();
		for(int i=0; i<hm.size(); i++)
			out.add(null);
		//TODO: use faster sorting algorithm
		Iterator<String> it = hm.keySet().iterator();
		while(it.hasNext()) {
			String s = it.next();
			Iterator<String> it2 = hm.keySet().iterator();
			int under = 0;
			while(it2.hasNext()) {
				String s2 = it2.next();
				if(s.equals(s2))
					continue;
				else if(hm.get(s) > hm.get(s2))
					under++;
			}
			while(out.get(under)!=null)
				under++;
			out.set(under, s);
		}
		return out;
	}
	
	public static HashMap<String, Integer> readDictionary(String inFile) throws Exception {
		BufferedReader brIn = new BufferedReader(new FileReader(inFile));
		String sLine;
		HashMap<String, Integer> dict = new HashMap<String, Integer>();
		while((sLine = brIn.readLine())!=null) {
			String [] fields = sLine.split("\t");
			dict.put(fields[0], Integer.parseInt(fields[1]));
		}
		return dict;
	}

	public static TreeSet<String> readTreeSet(String inFile) throws Exception {
		BufferedReader brIn = new BufferedReader(new FileReader(inFile));
		String sLine;
		TreeSet<String> out = new TreeSet<String>();
		while((sLine = brIn.readLine())!=null) {
			out.add(sLine);
		}
		return out;
	}

	
	public static void readArray(double[][] aO, BufferedReader brIn) throws
    Exception {
		for(int i=0; i < aO.length; i++) {
		    readArray(aO[i], brIn);
		}
	}
		public static void readArray(double[] aO, BufferedReader brIn) throws
		    Exception {
		for (int i = 0; i < aO.length; i++) {
		    String sLine = brIn.readLine();
		    aO[i] = Double.parseDouble(sLine);
		}
	}
		
	public static void writeArray(double[][] aO, BufferedWriter bwOut) throws
		    Exception {
		for (int i = 0; i < aO.length; i++) {
		    writeArray(aO[i], bwOut);
		}
	}
	public static void writeArray(double[] aO, BufferedWriter bwOut) throws
		    Exception {
		for (int i = 0; i < aO.length; i++)
		    bwOut.write(aO[i] + "\r\n");
	}
	
	public void renormalize() {
		for(int i=0; i<_pZgivenW.length; i++) {
			double sum = 0.0;
			for(int j=0; j<_pZgivenW[i].length; j++)
				sum += _pZgivenW[i][j];
			for(int j=0; j<_pZgivenW[i].length; j++)
				_pZgivenW[i][j] /= sum;

		}
		for(int i=0; i<_pCgivenZ.length; i++) {
			double sum = 0.0;
			for(int j=0; j<_pCgivenZ[i].length; j++)
				sum += _pCgivenZ[i][j];
			for(int j=0; j<_pCgivenZ[i].length; j++)
				_pCgivenZ[i][j] /= sum;
		}
	}
	
	//sets parameters randomly and re-normalizes
	public void randomInit() {
		for(int i=0; i<_pZgivenW.length; i++) 
			for(int j=0; j<_pZgivenW[i].length; j++)
				_pZgivenW[i][j] = Math.random();
		for(int i=0; i<_pCgivenZ.length; i++) 
			for(int j=0; j<_pCgivenZ[i].length; j++)
				_pCgivenZ[i][j] = Math.random();
		renormalize();
	}
	
	public static TextMasterEM readEMModel(String inFile) throws Exception {
		BufferedReader brIn = new BufferedReader(new FileReader(inFile));
		String sLine;
		sLine = brIn.readLine();
		int topics = Integer.parseInt(sLine);
		TextMasterEM outTM = new TextMasterEM(topics);
		readArray(outTM._pZgivenW, brIn);
		readArray(outTM._pCgivenZ, brIn);
		brIn.close();
		return outTM;
	}
	
	public void writeEMModel(String outFile, TextMasterEM outTM) throws Exception {
		BufferedWriter bwOut = new BufferedWriter(new FileWriter(outFile));
		String sLine;
		
		bwOut.write(outTM.TOPICS + "\r\n");
		writeArray(_pZgivenW, bwOut);
		writeArray(_pCgivenZ, bwOut);
		bwOut.close();
	}
	
	//HACK:hardcoded data length for this assignment
	public static int[][] readData(String inFile) throws Exception {
		BufferedReader brIn = new BufferedReader(new FileReader(inFile));
		String sLine;
		int [][] data = new int[121568][2];
		for(int i=0; i<data.length;i++) {
			sLine = brIn.readLine();
			String [] fields = sLine.split(" ");
			for(int j=0; j<2;j++)
				data[i][j] = Integer.parseInt(fields[j]);
		}
		brIn.close();
		return data;
	}
	
	public void trainEM(String dataFile, int ITERATIONS, String outputModelFile) throws Exception {
		int [][] data = readData(dataFile);
		randomInit();
		
		double [][] expectedWtoZ = new double[WORDS][TOPICS];
		double [][] expectedZtoC = new double[TOPICS][CONTEXTS];
		
		for(int iter=0; iter<ITERATIONS;iter++) {
			double logLikelihood = 0.0;
			//sum up expected sufficient statistics
			for(int i=0; i<data.length;i++) {
				int wordID = data[i][0];
				int contextID = data[i][1];
				
				//estimate P(Z | wordID, contextID):
				double [] pZCGivenW = new double[TOPICS];
				double sum = 0.0;
				

				//below, you must do three things:
				//1. increment logLikelihood 
				//  with the current conditional likelihood P(contextID | wordID) for this
				// example given the model
				//2. increment for each z expectedWtoZ with the expected number of observations of
				//  (wordID and Z) due to this example, and likewise 
				//3. increment for each z expectedZtoC with the expected number of observations of
				//  (Z and contextID) due to this example

				//NOTE: you can use the code I wrote below as a starting point, or not
				//The correct algorithm can be obtained by uncommenting and completing each 
				//of the four lines that start with //SUGGESTION:
				//The only variables you need are those local to this function, and the
				//class variables representing the parameters of the model ( _pZgivenW
				// and _pCgivenZ.

				//BEGIN code to be changed for assignment
				//First, set pZCGivenW equal to the likelihood of generating Z and C given W, i.e. P(C, Z | W).
				for(int z=0;z<TOPICS;z++) {
					//SUGGESTION: pZCGivenW[z] = ...
					sum += pZCGivenW[z];
				}
				//now increment the running log likelihood counter by log P(C | W):
				//SUGGESTION: logLikelihood += ...

				//now we normalize the P(C, Z | W) terms to sum to one for each Z, giving the expected number of times each Z occurs:
				double [] pZGivenCW = new double[TOPICS];
				for(int z=0;z<TOPICS;z++)
					pZGivenCW[z] = pZCGivenW[z] / sum;

				//finally, add up this example's contribution to the expected statistics:
				for(int z=0;z<TOPICS;z++) {
					//SUGGESTION: expectedWtoZ[wordID][z] += ...
					//SUGGESTION: expectedZtoC[z][contextID] += ...
				}

				//END code to be changed for assignment
			}
			//estimate new parameters based on expected sufficient statistics:
			for(int i=0; i<_pZgivenW.length; i++) 
				for(int j=0; j<_pZgivenW[i].length; j++)
					_pZgivenW[i][j] = expectedWtoZ[i][j];
			for(int i=0; i<_pCgivenZ.length; i++) 
				for(int j=0; j<_pCgivenZ[i].length; j++)
					_pCgivenZ[i][j] = expectedZtoC[i][j];
			renormalize();
			if(iter % 100 == 0)
				System.out.println("Iteration " + iter + " ALL: " + (logLikelihood/(double)data.length) + " " + Arrays.toString(_pZgivenW[0]));
		}
		writeEMModel(outputModelFile, this);
	}
	
    private static double KL(double [] a, double [] b) throws Exception
    {
        if(a.length != b.length)
            throw new Exception("Arrays must have same length.");
        double sum = 0;
        for(int i=0; i<a.length;i++)
        {
        	if(a[i]==0 || b[i] ==0)
        		continue;
            sum += a[i] * Math.log(a[i] / b[i]);
        }
        return sum;
    }

	
	public static void testWordVectors(String modelFile, String wordDictionary, String correctExamplesFile,
			String [] seedExamples) throws Exception {
		TextMasterEM tm = readEMModel(modelFile);
		HashMap<String, Integer> wordDict = readDictionary(wordDictionary);
		
		//build prototype array
		double [] ptype = new double[tm.TOPICS];
//		for(int j=0; j<tm.TOPICS;j++) {
//			ptype[j] = 1.0;
//		}
		double sum = 0.0;
		for(int i=0; i<seedExamples.length;i++) {
			int id = wordDict.get(seedExamples[i]);
			sum = 0.0;
			for(int j=0; j<tm.TOPICS;j++) {
				ptype[j] += tm._pZgivenW[id][j];
//				ptype[j] = Math.min(tm._pZgivenW[id][j], ptype[j]);
				sum+=ptype[j];
			}
		}
		for(int j=0; j<tm.TOPICS;j++) {
			ptype[j] /= sum;
			System.out.print(ptype[j] + " ");
		}
		System.out.println();
		
		//score examples by distance from prototype array
		TreeSet<String> correctExamples = readTreeSet(correctExamplesFile);
		HashMap<String, Double> wordScores = new HashMap<String, Double>();
		
		Iterator<String> wordIt = wordDict.keySet().iterator();
		while(wordIt.hasNext()) {
			String word = wordIt.next();
			int id = wordDict.get(word);
			double score = KL(tm._pZgivenW[id], ptype)+KL(ptype, tm._pZgivenW[id]);
			if(Double.isNaN(score))
				score = 1000000.0;
			wordScores.put(word, score);
		}
		
		//compute average precision
		TreeSet<String> seedTS = new TreeSet<String>();
		for(int i=0; i<seedExamples.length; i++)
			seedTS.add(seedExamples[i]);

		ArrayList<String> sortedWords = sortKeysByVal(wordScores);
		System.out.println("Sorted List: ");
		wordIt = sortedWords.iterator();
		while(wordIt.hasNext()) {
			String word = wordIt.next();
			System.out.println(word + "\t" + wordScores.get(word));
		}
		
		wordIt = sortedWords.iterator();
		double sumPrecision = 0.0;
		double currentPrecision = 0.0;
		double examplesProcessed = 0.0;
		double totalPositive = 0.0;
		
		while(wordIt.hasNext()) {
			String word = wordIt.next();
			if(seedTS.contains(word)) continue;
			if(correctExamples.contains(word)) {
				currentPrecision = (currentPrecision * (examplesProcessed) + 1.0) / (examplesProcessed + 1.0);
				sumPrecision += currentPrecision;
				totalPositive++;
			}
			else
				currentPrecision = (currentPrecision * (examplesProcessed) + 0.0) / (examplesProcessed + 1.0);
			examplesProcessed++;
			if(totalPositive >= (correctExamples.size() - seedTS.size()) - 0.5)
				break;
		}
		System.out.println("Average Precision " + (sumPrecision/totalPositive) + " vs. baseline " + ((double)(correctExamples.size()-seedTS.size())/(double)(sortedWords.size()-seedTS.size())));
	}
	
	public static void main(String [] args) throws Exception {
		try {
		if(args[0].equalsIgnoreCase("train"))
		{
			TextMasterEM tm = new TextMasterEM(Integer.parseInt(args[4]));
			tm.trainEM(args[1], Integer.parseInt(args[2]), args[3]);			
			return;
		}
		if(args[0].equalsIgnoreCase("test")) {
			String [] seeds = new String[] {args[4], args[5], args[6]};
			testWordVectors(args[1], args[2], args[3], seeds);
			return;
		}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		System.out.println("Usage: ");
		System.out.println("TextMasterEM train <data file> <num iterations> <output model file> <num topics>");
		System.out.println("TextMasterEM test <model file> <word dictionary> <correct examples file> <seed1> <seed2> <seed3>");
		
		
	}
}
