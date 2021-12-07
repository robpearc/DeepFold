/*******************************************************************************************
**  
**  Functions for loading the protein sequence and sec struct information
**
**  Please report bugs and questions to robpearc@umich.edu
**
*******************************************************************************************/

#ifndef PARSESEQ_H
#define PARSESEQ_H

#include "CommonParameters.h"
#include "Operations.h"

class ParseSeq  
{
public:
    ParseSeq();
    int seqnum;
    char *seqdata;
    char *seqheader;
    bool loadseq(const char *seqfile);
    bool loadseq2(const char *seqfile);
    bool loadss2(const char *ss2file);
    bool genesse();
    bool genesse2();
    void modifyss2();
    int getlongestH();
    virtual ~ParseSeq();
    
public: 
    double *dispos,*disposcp;
    ssef  *ss2;
    int numss2;
    int numsse;
    sssegment *sse;//for psipred
};

ParseSeq::ParseSeq()
{
    seqnum=0;
    numss2=0;
    numsse=0;
    sse=NULL;
    dispos=NULL;
    disposcp=NULL;
    ss2=NULL;
    seqheader=NULL;
    seqdata=NULL;
}

ParseSeq::~ParseSeq()
{
    if(seqdata)
    {
        delete[]seqdata;
        seqdata=NULL;
    }
    if(seqheader)
    {
        delete[]seqheader;
        seqheader=NULL;
    }
    if(ss2)
    {
        delete[]ss2;
        ss2=NULL;
    }
    if(sse)
    {
        delete[]sse;
        sse=NULL;
    }
    if(dispos)
    {
        delete[]dispos;
        dispos=NULL;
    }
    if(disposcp)
    {
        delete[]disposcp;
        disposcp=NULL;
    }
}

bool ParseSeq::loadseq2( const char *seqfile )
{
    FILE *file;
    char tmpc;
    file= fopen(seqfile, "rt");
    if(file==NULL)
    {
        fprintf(stderr,"Error when loading sequence file %s\n",seqfile);
        return false;
    }
    seqnum=0;
    if(!seqdata) seqdata=new char[65525];
    while(!feof(file))
    {
        tmpc=fgetc(file);
        if(feof(file))
        {
            fclose(file);
            seqdata[seqnum]='\0';
            return true;
        }
        if (!((tmpc>='A' && tmpc<='Z') || tmpc=='\r' || tmpc=='\n'))
            seqdata[seqnum++]='A';
        else if(tmpc>='A' && tmpc<='Z') seqdata[seqnum++]=tmpc;
    }
    fclose(file);
    seqdata[seqnum]='\0';
    return true;
}

/* read fasta or plain text protein sequence */
bool ParseSeq::loadseq(const char *seqfile)
{
    FILE *file;
    char tmpc;
    file= fopen(seqfile, "rt");
    if(file==NULL)
    {
        fprintf(stderr,"Error when loading sequence file %s\n",seqfile);
        return false;
    }
    seqnum=0;
    if(!seqheader) seqheader=new char[65525];
    if(!seqdata) seqdata=new char[65525];
    fgets(seqheader, 65525, file);
    if(seqheader[0]!='>')
    {
        delete[]seqdata;
        seqdata=NULL;
        delete[]seqheader;
        seqheader=NULL;
        fclose(file);
        return loadseq2(seqfile);
    }
    while(!feof(file))
    {
        tmpc=fgetc(file);
        if(feof(file))
        {
            fclose(file);
            seqdata[seqnum]='\0';
            return true;
        }
        if (!((tmpc>='A' && tmpc<='Z') || tmpc=='\r' || tmpc=='\n'))
            seqdata[seqnum++]='A';
        else if(tmpc>='A' && tmpc<='Z')
            seqdata[seqnum++]=tmpc;
    }
    fclose(file);
    seqdata[seqnum]='\0';
    return true;
}

bool ParseSeq::loadss2(const char *ss2file)
{
    FILE *file2;
    file2=fopen(ss2file,"rt");
    if(file2==NULL)
    {
        printf("Error when load ss2 file %s\n",ss2file);
        return false;
    }
    if(ss2)
    {
        delete[]ss2;
        ss2=NULL;
    }
    if(!seqheader) seqheader=new char[65525];
    if(!seqdata) seqdata=new char[65525];
    int allocss2=100;
    numss2=0;
    ss2=new ssef[allocss2];
    int i;
    char tmpstring[255];
    fgets(tmpstring, 255, file2);
    fgets(tmpstring, 255, file2);
    while(!feof(file2))
    {
        fgets(tmpstring, 255, file2);
        sscanf(tmpstring+5,"%c %c %f %f %f",&ss2[numss2].res,&ss2[numss2].ss,
               &ss2[numss2].a[0],&ss2[numss2].a[1],&ss2[numss2].a[2]);
        numss2++;
        if(numss2>=allocss2)
        {
            allocss2*=2;
            ss2=(ssef *)realloc(ss2,allocss2*sizeof(ssef));
        }
    }
    fclose(file2);
    numss2--;
    if(numss2!=0) ss2=(ssef *)realloc(ss2,numss2*sizeof(ssef));
    seqnum=numss2;
    for(i=0;i<numss2;i++) seqdata[i]=ss2[i].res;
    seqdata[numss2]='\0';
    //know the proba of disconnect
    if(dispos)
    {
        delete[]dispos;
        dispos=NULL;
    }
    dispos=new double[numss2];
    for(i=0;i<numss2-1;i++) dispos[i]=ss2[i].a[0]+ss2[i+1].a[0];
    dispos[numss2-1]=0;
    if(disposcp)
    {
        delete[]disposcp;
        disposcp=NULL;
    }
    disposcp=new double[numss2];
    for(i=0;i<numss2;i++) disposcp[i]=dispos[i];
    return true;
}

void ParseSeq::modifyss2()
{
    if(!ss2 || numss2==0) return;
    int i,j,k;
    for(i=1;i<numss2-1;i++)
    {
        if(ss2[i].ss=='H' && ss2[i-1].ss=='E' && ss2[i+1].ss=='E' && ss2[i].a[0]<ss2[i].a[2])
            ss2[i].ss='E';
        else if(ss2[i].ss=='H' && ss2[i-1].ss=='E' && ss2[i+1].ss=='E')
            ss2[i].ss='C';
    }
    for(i=0;i<numss2;i++)
    {
        j=i;
        while(j<numss2 && ss2[j].ss==ss2[i].ss) j++;
        if((j-i==2 || j-i==1) && ss2[i].ss=='H')
            for(k=i;k<j;k++) ss2[k].ss='C';
        i=j-1;
    }
}

bool ParseSeq::genesse()
{
    if(!ss2) return false;
    int i,j;
    if(sse)
    {
        delete[]sse;
        sse=NULL;
    }
    sse=new sssegment[numss2];
    numsse=0;
    for(i=0;i<numss2;i++)
    {
        j=i;
        while(j<numss2 && ss2[j].ss==ss2[i].ss) j++;
        if(j-i>=100 && ss2[i].ss=='C')
        {
            j=i+50;
            sse[numsse].init=i;
            sse[numsse].term=j-1;
            sse[numsse].ss=ss2[i].ss;
            numsse++;
            i=j-1;
        }
        else if(j-i>=3)
        {
            sse[numsse].init=i;
            sse[numsse].term=j-1;
            sse[numsse].ss=ss2[i].ss;
            numsse++;
            i=j-1;
        }
        else if(j-i==2)//two
        {
            if(i==0)//head then use the next one
            {
                j++;
                sse[numsse].init=i;
                sse[numsse].term=j-1;
                sse[numsse].ss=ss2[i].ss;
                numsse++;
                i=j-1;
            }
            else if(j==numss2)//tail
            {
                if(sse[numsse-1].term-sse[numsse-1].init==2)//previous has only two then merge the remaining two
                {
                    sse[numsse-1].term=numss2-1;
                    i=j-1;
                }
                else//use the previous one
                {
                    i--;
                    sse[numsse-1].term--;
                    sse[numsse].init=i;
                    sse[numsse].term=j-1;
                    sse[numsse].ss=ss2[i+1].ss;
                    numsse++;
                    i=j-1;
                }
            }
            else//middle
            {
                //BasicFunc bf;
                double tmp1,tmp2;
                tmp1=maxinthree(ss2[i-1].a[0],ss2[i-1].a[1],ss2[i-1].a[2]);
                tmp2=maxinthree(ss2[j].a[0],ss2[j].a[1],ss2[j].a[2]);
                if(tmp1<tmp2 && sse[numsse-1].term-sse[numsse-1].init>2)//use previous one
                {
                    i--;
                    sse[numsse-1].term--;
                    sse[numsse].init=i;
                    sse[numsse].term=j-1;
                    sse[numsse].ss=ss2[i+1].ss;
                    numsse++;
                    i=j-1;
                }
                else//use later one
                {
                    j++;
                    sse[numsse].init=i;
                    sse[numsse].term=j-1;
                    sse[numsse].ss=ss2[i].ss;
                    numsse++;
                    i=j-1;
                }
            }
        }
        else if(j-i==1)//one
        {
            if(i==0)//begin the first seg
            {
                if(ss2[j].ss!=ss2[j+1].ss) j++;
                while(j<numss2-1 && ss2[j].ss==ss2[j+1].ss) j++;
                j++;
                sse[numsse].init=i;
                sse[numsse].term=j-1;
                sse[numsse].ss=ss2[j-1].ss;
                numsse++;
                i=j-1;
            }
            else//merge to the former
            {
                sse[numsse-1].term=j-1;
                i=j-1;
            }
        }
    }
    if(numsse!=0) sse=(sssegment *)realloc(sse,numsse*sizeof(sssegment));
    
    //know the proba of disconnect
    if(dispos)
    {
        delete[]dispos;
        dispos=NULL;
    }
    dispos=new double[numss2];
    for(i=0;i<numss2-1;i++) dispos[i]=ss2[i].a[0]+ss2[i+1].a[0];
    for(i=0;i<numsse-1;i++) dispos[sse[i].term]+=2.0;
    dispos[numss2-1]=0;
    if(disposcp)
    {
        delete[]disposcp;
        disposcp=NULL;
    }
    disposcp=new double[numss2];
    for(i=0;i<numss2;i++)
    {
        disposcp[i]=dispos[i];
    }
    return true;
}

bool ParseSeq::genesse2()
{
    if(!ss2) return false;
    int i,j;
    if(sse)
    {
        delete[]sse;
        sse=NULL;
    }
    sse=new sssegment[numss2];
    numsse=0;
    for(i=0;i<numss2;i++)
    {
        j=i;
        while(j<numss2 && ss2[j].ss==ss2[i].ss) j++;
        if(j-i>=100 && ss2[i].ss=='C')
        {
            j=i+50;
            sse[numsse].init=i;
            sse[numsse].term=j-1;
            sse[numsse].ss=ss2[i].ss;
            numsse++;
            i=j-1;
        }
        else if(j-i>=1)
        {
            sse[numsse].init=i;
            sse[numsse].term=j-1;
            sse[numsse].ss=ss2[i].ss;
            numsse++;
            i=j-1;
        }
    }
    if(numsse!=0) sse=(sssegment *)realloc(sse,numsse*sizeof(sssegment));
    return true;
}

int ParseSeq::getlongestH()
{
    int numlong=0;
    int ti;
    int i;
    for(i=0;i<numsse;i++)
    {
        if(sse[i].ss=='H')
        {
            ti=sse[i].term-sse[i].init+1;
            if(ti>numlong) numlong=ti;
        }
    }
    return numlong;
}

#endif 
