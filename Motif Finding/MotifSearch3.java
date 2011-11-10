import java.io.*;
import java.util.*;
import java.math.*;

public class MotifSearch3
{ 
   public static void main(String[] args) 
   { 
		int numSeq;
		String fileName=args[0];
		numSeq=getNumSeq(fileName);
		int numRuns=20;
		int i;//loop variable
		int x;//loop variable
		int y;//loop variable
		int z;//loop variable
		int j;//loop variable
		int motifLength;
		int bestScore=99999999;
		int randomSeq1=0;
		int randomSeq2=1;
		int tempInt;//temp variable for the sort
		String tempString;//temp variable for sort
		
		int[] tempArray;//temp array for sort
		tempArray= new int[numSeq];
		String[] allSeq;//Array of All Sequences
		allSeq = new String[numSeq];//Declare Space for Array... modify to be able to accept any size input
		int[] motifPositions;//Positions of Current Run of MotifFinding
		motifPositions= new int[numSeq];
		
		String[] ConsensusMotif;//scores different motifs from the 20 iterations
		ConsensusMotif= new String[numRuns];
		int[] motifScore;//array to score all the different scores from the 20 iterations
		motifScore= new int[numRuns];
		int[][] positionVectors;
		positionVectors= new int[numRuns][numSeq];
		
		String MotifA;//first motif being tested in the greedy part
		String MotifB;//second motif in the greedy part
		Random generator= new Random();//set up the random number generator
		
		//get the information from the fasta... which is set as the first argument
		numSeq=getNumSeq(fileName);
		allSeq=getFasta(fileName, numSeq);
		//Second argument is the motifLength
		motifLength=Integer.parseInt(args[1]);
		
		
		//Best Motif needs to be set to a default... 
		String bestMotifA=allSeq[0].substring(0, motifLength);
		String bestMotifB=allSeq[1].substring(0, motifLength);
		
		
		//Do GreedySearch numRuns Times
		for(y=0; y<numRuns; y++)
		{
			bestScore=9999999;//reset Best Score
			randomSeq1=generator.nextInt(numSeq);
			randomSeq2=generator.nextInt(numSeq);
			//computeRandomB
			//Greedy Algorithm
			for(i=0;i<=allSeq[randomSeq1].length()-motifLength;i++)//iterates across one Sequence
			{
				for(x=0;x<=allSeq[randomSeq2].length()-motifLength;x++)//iterates across the other sequence
				{
					MotifA=allSeq[randomSeq1].substring(i, i+motifLength);
					MotifB=allSeq[randomSeq2].substring(x, x+motifLength);
					if(Hamming(MotifA,MotifB) <bestScore)
					{
					
						bestScore=Hamming(allSeq[randomSeq1].substring(i, i+motifLength),allSeq[randomSeq2].substring(x, x+motifLength));
						motifPositions[randomSeq1]=i;
						motifPositions[randomSeq2]=x;
						bestMotifA=MotifA;
						bestMotifB=MotifB;
					}
				}
			}
			
			
		
		
			for(i=0;i<numSeq;i++)//iterates across different Sequences
			{
				motifPositions[i]=0;//reset the position before checking which is best
				bestScore=9999999;//reset best score for next looping
				for(x=0;x<=allSeq[i].length()-motifLength;x++)//iterates across the sequence
				{
					
					MotifA=allSeq[i].substring(x, x+motifLength);
					
					if(Hamming(MotifA,bestMotifA)<bestScore)
					{
						bestScore=Hamming(MotifA,bestMotifA);
						motifPositions[i]=x;
					}
					else if(Hamming(MotifA,bestMotifB)<bestScore)
					{
						bestScore=Hamming(MotifA,bestMotifB);
						motifPositions[i]=x;
					}
				
				}
			}
			//Store Consensus for Current Iteration
			ConsensusMotif[y]=getConsensus(numSeq, motifLength,motifPositions,allSeq);	
			//Store Total Score for current Iteration
			motifScore[y]=allScore(allSeq,motifPositions,ConsensusMotif[y]);
			//Store Position Vector
			for (z=0;z<numSeq;z++)//loop to set positions
			{
			positionVectors[y][z]=motifPositions[z];	
			}
		}
		//bubble sort
		for(i=0;i<motifScore.length;i++)
		{
			for(j=i;j<motifScore.length;j++)
			{
				if(motifScore[i]>motifScore[j])
				{
					//switch score
					tempInt=motifScore[i];
					motifScore[i]=motifScore[j];
					motifScore[j]=tempInt;
					
					//switch Motif
					tempString=ConsensusMotif[i];
					ConsensusMotif[i]=ConsensusMotif[j];
					ConsensusMotif[j]=tempString;
					
					//Switch Position Vectors
					tempArray=positionVectors[i];
					positionVectors[i]=positionVectors[j];
					positionVectors[j]=tempArray;
					
				}
			}
		}
		try
		{
			FileWriter outputFile=new FileWriter("output.txt");
			PrintWriter out= new PrintWriter (outputFile);
			
			
			for(y=0; y<numRuns; y++)//output loop
			{
				out.print("Motif:"+ConsensusMotif[y]+"\t");
				out.print("PositionVector:");
				for(i=0;i<numSeq;i++)
				{
					out.print((positionVectors[y][i]+1)+",");
				}
				out.println("\t"+"Score:"+(float)motifScore[y]/numSeq);
			}
			out.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}		
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
	
	public static String getConsensus(int numSeq, int motifLength, int[] consensusPosition, String[] SeqArray)
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
	
	/*
	Hamming returns the hamming distance between two strings of equal length
	*/
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
	
	/*
	All Score Method Returns the score for a single motif against all sequences
	*/
	public static int allScore(String[] allSeq, int[] posVector, String motif)
	{
	int score=0;
	int i=0;//loop variable
	for (i=0; i<allSeq.length; i++)//for every sequence
	{
		score+=Hamming(allSeq[i].substring(posVector[i] , posVector[i]+motif.length()), motif);
	}
	
	return score;
	}
	public static int getNumSeq(String fileName)
	{
		int numSeq=0;
		String strLine;
		
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
					numSeq++;
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
		return numSeq;
   }
	
}	