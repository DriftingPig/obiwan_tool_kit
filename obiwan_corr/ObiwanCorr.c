//gcc ObiwanCorr.c -o ObiwanCorr -lm -I$HOME/.conda/envs/myenv/include -L$HOME/.conda/envs/myenv/lib
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define PI 3.1415926
#define TotBins 27   
double Separation(double sinthi1,double costhi1,double sinphi1, double cosphi1, double sinthi2, double costhi2, double sinphi2, double cosphi2){
    return costhi1*costhi2*(cosphi1*cosphi2+sinphi1*sinphi2)+sinthi1*sinthi2;
}

int getlinenum(char* f0){
    int ch, number_of_lines = 0;
    FILE* myfile;
    myfile = fopen(f0, "r");
    
    do 
    {
        ch = fgetc(myfile);
        if(ch == '\n')
        number_of_lines++;
    } while (ch != EOF);

    if(ch == '\n' && number_of_lines != 0) 
        number_of_lines--;

    fclose(myfile);

    return number_of_lines;
}

double** reading(char* f1)
{

    FILE *file;
    double** data;
    double test1,test2,test3,test4;
    int linenum=0;
    int i,j;
    linenum = getlinenum(f1);
    file = fopen(f1, "r");
    data = (double **) malloc(linenum*(sizeof(double *)));
    for (i=0; i<linenum; i++)
    {
        data[i] = (double*) malloc(5*(sizeof(double)));
    }
    
    for (i = 0; i < linenum; i++)
        for (j=0; j<5; j++)
    {
        fscanf(file, "%lf", &data[i][j]);
    }
    fclose(file);
    return data;
}

double* BinAllocator(int bins){
    double MaxAngle = log10(5.*PI/180.);
    double MinAngle = log10(0.01*PI/180.);
    double AngleInteval = (MaxAngle-MinAngle)/(float)bins;
    int i;
    double* binns;
    binns = (double *)malloc((bins+1)*sizeof(double));
    for(i=0;i<bins+1;i++){
    binns[i] = cos(pow(10,i*AngleInteval+MinAngle));
    }
    return binns;
}

int BinNum(double AngD, double Binnings[],int Totbins){
    int ii;
    int count=0;
    int flag=0;
    for(ii=1;ii<Totbins+1;ii++){
    if(AngD>=Binnings[ii]){
        flag=1;
        break;}
        
    }
    if(AngD>=Binnings[0]){
    return -2;
    }
    if(flag){
    return ii;
    }
    else{
    return -1;
    }
}
//pairtype:DD 0,DR 1, RR 2
//./Obiwan_corr filename1 filename2 self_flag outputfilename mode pairtype
int main(int argc, char* argv[]){
    char* filename1 = argv[1];
    char* filename2 = argv[2];
    char* outputfilename = argv[4];
    int self_flag = atoi(argv[3]), pairtype = atoi(argv[6]), index1 = atoi(argv[7]),index2 = atoi(argv[8]),bin_flag=0;
    char* mode = argv[5];
    FILE *file1, *file2, *FileOut;
    double** dataset1;
    double** dataset2;
    int Bin_Number;
    double outsidebin=0,smallsidebin=0,totalpts=0;
    static double HistBins[TotBins];
    double* binning;
    double AngDis;
    int i,j,TotLine1,TotLine2;
    //read in all the data
    dataset1 = reading(filename1);
    dataset2 = reading(filename2);
    //set up the binning
    binning = BinAllocator(TotBins);
    //get the line number
    TotLine1 = getlinenum(filename1);
    TotLine2 = getlinenum(filename2);
    //calculate the angular distance, put a number in the corresponding bin
    if(!self_flag){
    for(i=0;i<TotLine1;i++)
        for(j=0;j<TotLine2;j++){
    AngDis = Separation(dataset1[i][0],dataset1[i][1],dataset1[i][2],dataset1[i][3],dataset2[j][0], dataset2[j][1], dataset2[j][2], dataset2[j][3]);
    totalpts+=dataset1[i][4]*dataset2[j][4];
    Bin_Number = BinNum(AngDis, binning, TotBins);
    if(Bin_Number==-1){outsidebin+=dataset1[i][4]*dataset2[j][4];}
    if(Bin_Number==-2){smallsidebin+=dataset1[i][4]*dataset2[j][4];}
    else if(Bin_Number<=TotBins){
        HistBins[Bin_Number-1]+=dataset1[i][4]*dataset2[j][4];}
    }
    }
    else{
    for(i=0;i<TotLine1;i++)
        for(j=i+1;j<TotLine2;j++){
    AngDis = Separation(dataset1[i][0],dataset1[i][1],dataset1[i][2],dataset1[i][3],dataset2[j][0], dataset2[j][1], dataset2[j][2], dataset2[j][3]);
    totalpts+=dataset1[i][4]*dataset2[j][4];
    Bin_Number = BinNum(AngDis, binning, TotBins);
    if(Bin_Number==-1){outsidebin+=dataset1[i][4]*dataset2[j][4];}
    if(Bin_Number==-2){smallsidebin+=dataset1[i][4]*dataset2[j][4];}
    else if(Bin_Number<=TotBins){
        HistBins[Bin_Number-1]+=dataset1[i][4]*dataset2[j][4];}
    }    
    }
    printf("%lf pairs outside the bin, %lf pairs too small\n",outsidebin,smallsidebin);
    FileOut = fopen(outputfilename, "w");
    for(i=0;i<TotBins;i++){
    fprintf(FileOut,"%d %lf\n",i,HistBins[i]);
    }
    fclose(FileOut);
    FileOut = fopen("/global/cscratch1/sd/huikong/obiwan_data/data_from_obiwantest/obiwan_test/obiwan_corr/data/BinHist/weighed_uniform/TotalPoints.txt",mode);
    if(pairtype==0){
    fprintf(FileOut,"%lf 0 0 %d %d\n",totalpts,index1,index2);}
    if(pairtype==1){
    fprintf(FileOut,"0 %lf 0 %d %d\n",totalpts,index1,index2);}
    if(pairtype==2){
    fprintf(FileOut,"0 0 %lf %d %d\n",totalpts,index1,index2);}
    fclose(FileOut);
    return(0);
}
