package kclust;
import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.lang.*;

public class k_means_clustering
{
  public static void main (String args[]) throws IOException
	{

	 Centroid cent = new Centroid();
	 int ClustNumber;
	 System.out.println(" Enter the number of clusters");
	 Scanner input = new Scanner(System.in);
	 ClustNumber=input.nextInt();
	 String [][] numbers = new String [6][2];
	 double Cordx[] =new double[6];
	 double Cordy[] =new double[6];
	 File file = new File("/Users/tulip_nandu/Documents/Algorithms/resources/sam.csv");
	 BufferedReader bufRdr = new BufferedReader(new FileReader(file));
	 String line = null;
	 int row = 0;
	 int col = 0;

	 //read each line of text file
	 while((line = bufRdr.readLine()) != null && row< 6 )
	 {
	 StringTokenizer st = new StringTokenizer(line,",");
	 while (st.hasMoreTokens())
	 {
	 //get next token and store it in the array
	 numbers[row][col] = st.nextToken();
	 col++;
	 }
	 col = 0;
	 row++;
	 }


	 for(row=0;row < 6;row++)
	 {
	 for(col=0; col<2;col++)
	 {
	 System.out.print(" " + numbers[row][col]);
	 }
	 System.out.println(" ");
	 }
	 for(row=0;row<6;row++)
	 {

	 Cordx[row]=Double.parseDouble(numbers[row][0]);
	 Cordy[row]=Double.parseDouble(numbers[row][1]);
	 //System.out.print(" " + Cordx[row]);
	 }

	 for(row=0;row<6;row++)
	 {
	 System.out.print(" " + Cordx[row]);
	 //System.out.print("\n " + Cordy[row]);
	 }
	 System.out.print(" \n");
	 for(row=0;row<6;row++)
	 {

	 System.out.print(" " + Cordy[row]);
	 }
	 cent.Grouping(Cordx,Cordy,ClustNumber);

	 }
	 }

