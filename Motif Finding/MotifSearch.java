import java.io.*;
import java.math.*;

public class MotifSearch 
{ 
   public static void main(String[] args) 
   { 
		int numSeq=5;
		int lengthSeq=10;
		int i;//loop variable
		int motifLength=5;
		String[] allSeq;//Array of All Sequences
		allSeq = new String[numSeq];//Declare Space for Array... modify to be able to accept any size input
		int[] motifPositions;//Positions of Current Run of MotifFinding
		motifPositions= new int[numSeq];
		String ConsensusMotif="ggagt";
		
		
		
		
		
		allSeq=getFasta(args[0], numSeq);
		
		motifPositions[0]=5;
		motifPositions[1]=0;
		motifPositions[2]=2;
		motifPositions[3]=1;
		motifPositions[4]=1;
		System.out.println("Score "+ allScore(allSeq, motifPositions, ConsensusMotif ));
		//for(i=0;i<numSeq;i++)//this is for testing... sets all positions to 1
		//{motifPositions[i]=1;}
		
		ConsensusMotif=getConsensus(numSeq,lengthSeq, motifLength,motifPositions,allSeq);
		System.out.println("Consensus is "+ConsensusMotif);
		
   }
   /*
   METHOD: getFasta
   ARGUMENTS: File Name to be read from, the number of sequences to be read
   RETURNS: An array of strings, each sequence in a different element of the array
   
   Sets up a buffered reader to read from the fasta file.  It reads the fasta as a string
   one line at a time.  It will throw out the line if it detects a > character indicating
   the string is not a sequence. It will then iterate thru each string and put it into the array
   using toLowerCase to make sure all characters can be compared to one another. Streams are then closed.
   */
   
   public static String[] getFasta(String fileName, int numSeq)
   {
   
   int x=0;//for the loop
   String strLine;
   String[] SeqArray;
   SeqArray= new String[numSeq];
    
   try
		{
			
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			
			while ((strLine = br.readLine()) != null)
			{
				//Check to see if the Line is Sequence or 
				if((strLine.charAt(0)!=('>')))//all sequence lines start with >
				{
								
						SeqArray[x]=strLine.toLowerCase();
					
					x++;//moves onto the next sequence
				}
			}
			
			//close Streams
			br.close();
			fstream.close();
			in.close();
		}
		catch (Exception e)
		{
			System.err.println("Error: " + e.getMessage());
		}
		return SeqArray;
	}
	/*
	METHOD: getConsensus
	ARGUMENTS: Number of Sequences, Length of sequences, length of Motif, a array of consensus Positions, an array of all sequences
	RETURNS: a consensus string 
	*/
	
	public static String getConsensus(int numSeq, int lengthSeq, int motifLength, int[] consensusPosition, String[] SeqArray)
	{
		String ConsensusString="";
		int i;//iterates thru sequences
		int x;//iterates thru characters
		int aCount=0;
		int tCount=0;
		int gCount=0;
		int cCount=0;
		int currentPostion;
		char currentChar='a';//default character
		
		
		
		for (x=0; x<motifLength; x++)
		{
			aCount=0;//we need to reset all counts before moving
			tCount=0;//to the next set of counts
			gCount=0;
			cCount=0;
			for(i=0; i<numSeq;i++)//go thru all sequences 
			{
				
				if(SeqArray[i].charAt(consensusPosition[i]+x)=='a')
				{aCount++;
					if (aCount>Math.max(tCount, Math.max(gCount,cCount)))
					{currentChar='a';}
				}
				if(SeqArray[i].charAt(consensusPosition[i]+x)=='t')
				{tCount++;
					if (tCount>Math.max(aCount,Math.max(gCount,cCount)))
					{currentChar='t';}
				}
				if(SeqArray[i].charAt(consensusPosition[i]+x)=='g')
				{gCount++;
					if (gCount>Math.max(tCount,Math.max(aCount,cCount)))
					{currentChar='g';}
				}
				if(SeqArray[i].charAt(consensusPosition[i]+x)=='c')
				{cCount++;
					if (cCount>Math.max(tCount,Math.max(gCount,aCount)))
					{currentChar='c';}
				}	
			}//Iterates Thru Sequences Loop Ends
			ConsensusString +=currentChar;
		}
		return ConsensusString;
	}
	public static int Hamming(String compOne, String compTwo)
	{
		int counter = 0;

		for (int i = 0; i < compOne.length(); i++)
		{
			if (compOne.charAt(i) != compTwo.charAt(i))
			{
				counter++;
			}
		}

		return counter;
	}
	
	
	public static int allScore(String[] allSeq, int[] posVector, String motif)
	{
	int score=0;
	int i=0;//loop variable
	System.out.println("HELLOOOO");
	for (i=0; i<allSeq.length; i++)//for every sequence
	{
		System.out.println(allSeq[i].substring(posVector[i], posVector[i]+motif.length()));
		score+=Hamming(allSeq[i].substring(posVector[i] , posVector[i]+motif.length()), motif);
	}
	
	return score;
	}
}	