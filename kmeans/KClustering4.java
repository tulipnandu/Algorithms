import java.util.Random;
import java.io.*;
import java.lang.Integer;

public class KClustering4
{
	public static void main(String[] args)
	{
		//Puts the FASTA into the array allSeq
		int numSeq;
		String fileName=args[0];//file name is first argument	
		numSeq=getNumSeq(fileName);
		String[] allSeq;//array of all sequences
		allSeq= new String[numSeq];
		allSeq=getFasta(fileName, numSeq);
		
		//set K as numClust.  if this value is not an int the program will not run
		int numClust;
		numClust=Integer.parseInt(args[1]);
		
		//Other Variables
		int[] seqNumbers;//array of the initial sequence numbers
		seqNumbers= new int[numClust];//can be random or set in file
		
		String[] seqNumString;//stores the numbers for the clusters
		seqNumString= new String[numClust];
		
		int x;//loop variable
		int y;//loop variable
		int z;//loop variable
		int q;//loop variable
		int i;//loop variable
		int score=0;//score for current score in consensus
		int tempVal;//for the bubble sort
		
		int lowScore[];//score for current best score in consensus
		lowScore= new int[numClust];
		
		
		int[] whatClust;
		whatClust= new int[numSeq]; //stores what values are in what cluster
		
		boolean hasChanged=false;//a variable that stores if a change has been made to the clusters
		
		int numIteration=0;
		String alphabet="atgc";//we could change this to protein sequence if we wanted to.
		
		int[] baseCount;
		baseCount= new int[alphabet.length()];
		
		char currentChar;
		int haveSeedFile=0;
		
		
		
		//find seed file if it exists
		String seedFile="";
		
		try
		{
			
			seedFile=args[2];
				
			
			seqNumString=getFasta(seedFile, numClust);//this method will give the numbers in the file without having to rewrite anything
			
			for(x=0; x<=numClust;x++)
			{
				
				seqNumbers[x]=Integer.valueOf(seqNumString[x]);//we do need to parse the string to an integer.
				
				haveSeedFile=1;
			}
			
			
		}
		catch(ArrayIndexOutOfBoundsException e)
		{
			
			if(haveSeedFile==0)
			{
				Random generator = new Random();
			
				for(x=0; x<numClust;x++)
				{
				seqNumbers[x]=generator.nextInt(numSeq);
				}
				
			}
		}

	
		do
		{
			hasChanged=false;
			for(x=0; x<numSeq; x++)
			{
				
				for(y=0;y<numClust;y++)
				{
					
					if(Hamming(allSeq[x],allSeq[seqNumbers[y]])<Hamming(allSeq[x],allSeq[seqNumbers[whatClust[x]]]))
					{
						whatClust[x]=y;
						
						hasChanged=true;// if this is not set the program will loop again
					}
				}
			}
			
			//computing new consensus
			for(x=0; x<numClust; x++)
			{
				
				lowScore[x]=99999999;
				
				for (y=0; y<allSeq.length; y++)//across each sequence
				{
					score=99999999;//score set high to avoid bad comparrison.
					if(whatClust[y]==x)
					{
						score=0;
						
						{
							for(z=0;z<allSeq.length; z++)
							{
								if(whatClust[z]==x)
								{
								score+=Hamming(allSeq[y], allSeq[z] );
								
								}
							}
						}	
					}	
					if(score<lowScore[x])
					{
						
						lowScore[x]=score;
						seqNumbers[x]=y;
						
						
					}	
					
				}
			}
		numIteration++;
		
		}//end k-means loop
		while(hasChanged==true);
		
		
		int[] clusID;
		clusID= new int[numClust];
		for(y=0; y<numClust;y++)
		{
			clusID[y]=y;
		}	
		//sort output
		for(y=0; y<numClust;y++)
		{
			
			for(i=y; i<numClust;i++)
			{
				if(seqNumbers[y]>seqNumbers[i])
				{
					
					//swap seed ID values
					tempVal=seqNumbers[y];
					seqNumbers[y]=seqNumbers[i];
					seqNumbers[i]=tempVal;
					
					
					//Swap Cluster numbers
					tempVal=clusID[y];
					clusID[y]=clusID[i];
					clusID[i]=tempVal;
									
					
					//Swap Score
					tempVal=lowScore[y];
					lowScore[y]=lowScore[i];
					lowScore[i]=tempVal;
				}
			}
		}
		
		//output results to file
		try
		{
			FileWriter outputFile=new FileWriter("output.txt");
			PrintWriter out= new PrintWriter (outputFile);
			
			
			for(y=0; y<numClust; y++)//output loop for each cluster
			{
				out.print(seqNumbers[y]+"\t");
				
				for(i=0;i<numSeq;i++)
				{
					if(clusID[whatClust[i]]==y & i!=seqNumbers[y])
					{
					out.print(i+",");
					}
				}
				out.println("\t"+"Score:"+lowScore[y]);
			}
			out.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
		System.out.println("Output in file: output.txt");
	}//end main
	
	/*
   METHOD: getFasta
   ARGUMENTS: File Name to be read from, the number of sequences to be read
   RETURNS: An array of strings, each sequence in a different element of the array
   
   Sets up a buffered reader to read from the fasta file.  It reads the fasta as a string
   one line at a time.  It will throw out the line if it detects a > character indicating
   the string is not a sequence. It will then iterate thru each string and put it into the array
   using toLowerCase to make sure all characters can be compared to one another. Streams are then closed.
   
   This is the same method used from the first assignment
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
	Hamming returns the hamming distance between two strings of equal length
	Reused from assignment 1
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
	Returns number of sequences in a fasta
	Reused from first program
	*/
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