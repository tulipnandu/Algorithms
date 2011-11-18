package kclust;
import java.text.DecimalFormat;


public class Centroid
	{

	// finds the mean of the dataset for grouping
	public void Grouping(double[] Cordx, double[] Cordy, int clustNumber)
	{
	int clusterNumber = clustNumber;
	double[] ClustCordX = new double[clustNumber];
	double[] ClustCordY = new double[clustNumber];
	this.getMeansetCentroid(Cordx, Cordy, clustNumber);
	DecimalFormat dec = new DecimalFormat("0.00");
	for(int i = 0;i<Cordx.length;i++)
	{
	String result1 = dec.format(Cordx[i]);
	String result2 = dec.format(Cordy[i]);
	System.out.println("\n Cords are ( " + result1 + " , " + result2 + ")");
	}
	//setting random datasets as centroids
	for(int i = 0; i<clustNumber;i++)
	{
	ClustCordX[i] = Cordx[i];
	ClustCordY[i] = Cordy[i];
	}


	this.groupCordtoCluster(Cordx,Cordy,ClustCordX,ClustCordY);
	}

	public void groupCordtoCluster(double[] Cordx, double[] Cordy, double[] ClustCordX, double[] ClustCordY)
	{
	double temp ;
	int size = Cordx.length;
	int clustsize = ClustCordX.length;
	int clusterComparison = clustsize;
	int[] grouping = new int[size - clustsize];
	double[] ClustgroupX = new double[size - clustsize];
	double[] ClustgroupY = new double[size - clustsize];
	int tempint = -1;

	//grouping the dataset to respective clusters by comparing the distance
	for(int i = clusterComparison; i < size;i++)
	{
	temp = 0;
	for(int j = 0;j<clustsize;j++)
	{
	if (j == 0)
	tempint++;
	if(temp == 0)
	{
	temp = Math.sqrt(Math.pow((Cordx[i]-ClustCordX[j]),2) + Math.pow((Cordy[i]-ClustCordY[j]),2));
	grouping[tempint] = j;
	ClustgroupX[tempint] = Cordx[i];
	ClustgroupY[tempint] = Cordy[i];
	}
	else if (temp > Math.sqrt(Math.pow((Cordx[i]-ClustCordX[j]),2) + Math.pow((Cordy[i]-ClustCordY[j]),2)))
	{
	temp = Math.sqrt(Math.pow((Cordx[i]-ClustCordX[j]),2) + Math.pow((Cordy[i]-ClustCordY[j]),2));
	grouping[tempint] = j;
	ClustgroupX[tempint] = Cordx[i];
	ClustgroupY[tempint] = Cordy[i];
	}
	}

	}
	DecimalFormat dec = new DecimalFormat("0.00");
	String result1, result2, result3, result4;
	for(int i = 0; i<grouping.length;i++)
	{
	//for(int i = 1; i< clustNumber;i++) {
	System.out.println("------------------------");
	System.out.println("Clusters for group " + grouping[i]);
	result1 = dec.format(Cordx[grouping[i]]);
	result2 = dec.format(Cordy[grouping[i]]);
	result3 = dec.format(ClustgroupX[i]);
	result4 = dec.format(ClustgroupY[i]);
	System.out.println("Cordinates are (" + result1 + " , " + result2 + ")");
	System.out.println("------------------------");
	System.out.println("Clusters for group " + grouping[i]);
	System.out.println("Cordinates are (" + result3 + " , " + result4 + ")");

	}
	}

	public void getMeansetCentroid(double[] Cordx, double[] Cordy, int ClustNumber)
	{
	double xCord, yCord;
	double MAX=0,distance, tempd1, tempd2;
	double[] Distances = new double[Cordx.length];
	int reference, i, j, temp1, temp2, point, length;
	reference = i = j = temp1 = temp2 = point =0;
	int[] centroids;


	for(j = 1; j < Cordx.length;j++)
	{
	Distances[j-1] = Math.sqrt(Math.pow((Cordx[j]-Cordx[0]),2) + Math.pow((Cordy[j]-Cordy[0]),2));
	}
	for(i=0;i<Cordx.length-1;i++)
	{
	for(j=0;j<Cordx.length-1-i;j++)
	{
	if(Distances[j+1] < Distances[j])
	{
	distance = Distances[j];
	tempd1 = Cordx[j];
	tempd2 = Cordy[j];
	Distances[j] = Distances[j+1];
	Cordx[j] = Cordx[j+1];
	Cordy[j] = Cordy[j+1];
	Distances[j+1] = distance;
	Cordx[j+1] = tempd1;
	Cordy[j+1] = tempd2;
	}
	}
	}
	//recalculation of centroids
	point = Cordx.length;
	do
	{
	if(Cordx.length % ClustNumber != 0)
	point--;
	}while(point % ClustNumber != 0);
	length = point/ClustNumber;
	for(i=0;i<Cordx.length;i=length+i)
	{
	if((i+length-1) > point)
	break;
	tempd1 = Cordx[i];
	tempd2 = Cordy[i];
	Cordx[i] = Cordx[i+length-1];
	Cordy[i] = Cordy[i+length-1];
	Cordx[i+length-1] = tempd1;
	Cordy[i+length-1] = tempd2;
	}

	}


	}

