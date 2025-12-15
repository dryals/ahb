 #include<iostream>
 #include <string>
 #include <fstream>
 #include <cmath>
 #include <sstream>
 #include <cstdlib>
 #include <vector>
 #include <algorithm>
 using namespace std;

 int main(int argc, char** argv){
 { 
  if (argc < 5)
   {
   cerr << "Inputs: obsGenoFileName impGenoFileName outPutFileName numberOfcolumns" << endl;
   exit(1);
   }
   string obsGenoFileName = argv[1];
   string impGenoFileName = argv[2];
   string outPutFileName  = argv[3];
   const int numberOfColumns = atoi(argv[4]);
     
   ifstream inFile1(obsGenoFileName.c_str()); 
   ifstream inFile2(impGenoFileName.c_str()); 
   ofstream outFile(outPutFileName.c_str()); 
   
   vector <int> a (numberOfColumns);
   vector <int> b (numberOfColumns);

// Declare r valiables
   
vector <float> sx(numberOfColumns),sy(numberOfColumns),sxy(numberOfColumns),sxx(numberOfColumns),syy(numberOfColumns),mx(numberOfColumns),my(numberOfColumns),sdx(numberOfColumns),sdy(numberOfColumns),cxy(numberOfColumns),r(numberOfColumns),vx(numberOfColumns),vy(numberOfColumns),fr(numberOfColumns),accuracy(numberOfColumns),conrd(numberOfColumns);
// initialize a vector
   for (int i =0; i< numberOfColumns; i++)
    {
      fr[i]=0;
      a[i]=0;
      b[i]=0;
      conrd[i]=0;
      accuracy[i]=0;
      sx[i]=0;
      sy[i]=0;
      sxy[i]=0;
      sxx[i]=0;
      syy[i]=0;
      mx[i]=0;
      my[i]=0;
      vx[i]=0;
      vy[i]=0;
      sdx[i]=0;
      sdy[i]=0;
      cxy[i]=0;
      r[i]=0;
    }
   int val1,val2;
   int count =0;  
	int status = 1;
	while (status == 1)
	 {
		for (int i = 0; i < numberOfColumns; i++)
		{
			inFile1 >> val1;
	      inFile2 >> val2;
   		if (inFile1.good() != 1)
			{
				status = 0;
				break;
			}
 if(status==1) count++;
 if (val1!=5){
     sx[i]+=val1;
     sy[i]+=val2;
     sxx[i]+=val1*val1;
     syy[i]+=val2*val2;
     sxy[i]+=val1*val2;
     fr[i]++;
           }
      if(val1==val2){
      conrd[i]++;   
         }
      }
   }
 for (int i = 0; i < numberOfColumns; i++){
     mx[i]=sx[i]/fr[i]; 
     my[i]=sy[i]/fr[i];
     vx[i]=(sxx[i]/fr[i])-(mx[i])*(mx[i]);
     vy[i]=(syy[i]/fr[i])-(my[i]*my[i]);
     sdx[i]=sqrt(vx[i]);
     sdy[i]=sqrt(vy[i]);
     cxy[i]=(sxy[i]/fr[i])-(mx[i]* my[i]);
     r[i]= cxy[i]/(sdx[i]*sdy[i]);
     accuracy[i]=conrd[i]/fr[i];
    }

 outFile << "SNPID\taccuracy\tcorrelation" << endl;
 for (int i = 0; i < numberOfColumns; i++) 
	{

 outFile << i+1 << '\t' << accuracy[i] << '\t' << r[i] << endl;

}
	inFile1.close();
	inFile2.close();
	outFile.close();
}
  return 0;
}

