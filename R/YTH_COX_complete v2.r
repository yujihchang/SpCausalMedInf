library("survival")


YTH_COX_diff=function(time,A,X1,X2,D,X,expor=c(za1,za2,zb1,zb2),center="TRUE",setting=NULL,EFFECT,saving="FALSE"){
#time=P.time;X=matrix(Za);expor=c(1,2,1,2);EFFECT=c("DE","IE");center="TRUE";setting=NULL;saving="FALSE"


#t=tt; A=dat[,1]*0; X1=dat[,1]; X2=dat[,3]; D=dat[,4];X=as.matrix(x);expor=c(0,1,0,1);EFFECT=c("IE","DE");center="TRUE";setting=NULL;saving="FALSE"



DE.out=IE.out=DE.sd.out=IE.sd.out=NULL
za1=expor[1];za2=expor[2];zb1=expor[3];zb2=expor[4]
n=length(X2)
DL1X=DL0X=rep(0,n)
dLan1X=dLan0X=0
#########################################################
SX2=sort(X2,index.return=TRUE)
X1=X1[SX2$ix]
A=A[SX2$ix]
X=as.matrix(X[SX2$ix,])
Dperp=D=D[SX2$ix]
X2=SX2$x
########################################################
#center="TRUE";X=as.matrix(XM);expor=c(1,2,1,2);EFFECT=c("DE","IE");t=time;expor=c(1,2,1,2)


X1M=matrix(X1,n,n)
X2M=matrix(X2,n,n)
N1t=(X1M<t(X2M))&(X1M<X2M)
p=NCOL(X)
dN21=diag(N1t)
if(saving=="TRUE"){gc()}
#########################################################

Y21=X2[dN21==1]
n21=length(Y21)

Y11=X1[dN21==1]
A21=A[dN21==1]
X21=matrix(X[dN21==1,],sum(dN21==1),p)
D21=D[dN21==1]

AY11=ifelse(Y11>= A21,Y11,A21)
B21=coxph(Surv(AY11,Y21,D21)~X21)
MdifY=matrix(Y21,n21,n21)
MadjA21=matrix(AY11,n21,n21)
diff_Risk_n1=( (MadjA21<=t(MdifY)&t(MdifY)<=MdifY) )

SXB21=exp(apply(t(X21)*B21$coef,2,sum))
diff_dR1=D21/apply(diff_Risk_n1*SXB21,2,sum)
###################################################################################
DL1X_za2=DL1X
if(center=="FALSE") sb1_za2=sum(B21$coef*c(za2,rep(0,(p-1))) )
if(center=="TRUE") sb1_za2=sum(B21$coef*c(za2,apply(as.matrix(X[,-1]),2,mean)))
if(center=="SET")  sb1_za2=sum(B21$coef*c(za2,setting))

dLan1X_za2=diff_dR1*exp(sb1_za2)
DL1X_za2[dN21==1]=dLan1X_za2
#SSS=sort(X2,index.return=TRUE);plot(SSS$x,cumsum(DL1X[SSS$ix]),type="s")
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##
if(sum(EFFECT%in%"DE")==1){
DL1X_za1=DL1X
if(center=="FALSE") sb1_za1=sum(B21$coef*c(za1,rep(0,(p-1))) )
if(center=="TRUE") sb1_za1=sum(B21$coef*c(za1,apply(as.matrix(X[,-1]),2,mean)))
if(center=="SET")  sb1_za1=sum(B21$coef*c(za1,setting))
dLan1X_za1=diff_dR1*exp(sb1_za1)
DL1X_za1[dN21==1]=dLan1X_za1
#SSS=sort(X2,index.return=TRUE);plot(SSS$x,cumsum(DL1X[SSS$ix]),type="s")
 }
#################################################################################################
if(saving=="TRUE"){gc()}
DL0X_za2=DL0X_za1=DL0X
Y20=X2
Y20[dN21==1]=X1[dN21==1]
##

#
B20=NULL;D20=NULL;diff_dR0=NULL;X20=A20=NULL
##
if(length(Y20)>0){ #預防沒有 workL
#dLan0X=rep(0,length(Y20))
D20=D
D20[dN21==1]=0


workL=which(A<Y20)
Y20=Y20[workL]
Y10=X1[workL]
A20=A[workL]
if( length(Y20)==1 ) X20=matrix(X[workL,],length(length(Y20)),p)
if( length(Y20)> 1 ) X20=as.matrix(X[workL,])
D20=D20[workL]


Y20=sort(Y20,index.return=TRUE)
D20=D20[Y20$ix]
A20=A20[Y20$ix]
X20=as.matrix(X20[Y20$ix,])
workL0=workL[Y20$ix]
Y20=Y20$x

B20=coxph(Surv(A20,Y20,D20)~X20)

n20=length(Y20)
MY20=matrix(Y20,n20,n20)
MA20=matrix(A20,n20,n20)
diff_Risk_n0=( (MA20<=t(MY20)&t(MY20)<=MY20) )
SXB20=exp(apply(t(X20)*B20$coef,2,sum))
diff_dR0=D20/apply(diff_Risk_n0*SXB20,2,sum)


#################################################################################################
if(saving=="TRUE"){gc()}

if(center=="FALSE")sb0_za2=sum(B20$coef*c(za2,rep(0,(p-1))) )
if(center=="TRUE") sb0_za2=sum(B20$coef*c(za2,apply(as.matrix(X[,-1]),2,mean)))
if(center=="SET") sb0_za2=sum(B20$coef*c(za2,setting))


dLan0X_za2=diff_dR0*exp(sb0_za2)
dLan0X_za2=ifelse(is.na(dLan0X_za2),0,dLan0X_za2)
DL0X_za2[workL0]=dLan0X_za2


#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##
if(sum(EFFECT%in%"DE")==1){ 
DL0X_za1=DL0X
if(center=="FALSE")sb0_za1=sum(B20$coef*c(za1,rep(0,(p-1))) )
if(center=="TRUE") sb0_za1=sum(B20$coef*c(za1,apply(as.matrix(X[,-1]),2,mean)))
if(center=="SET") sb0_za1=sum(B20$coef*c(za1,setting))
dLan0X_za1=diff_dR0*exp(sb0_za1)
dLan0X_za1=ifelse(is.na(dLan0X_za1),0,dLan0X_za1)
DL0X_za1[workL0]=dLan0X_za1
}
#-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##-----##
   }  #  for " if "預防沒有 workL
########################################################################################################
if(saving=="TRUE"){gc()}

AYm=ifelse(Y11>A21,Y11,A21)
adjA=A
#adjA[dN21==1]=AYm

AM=matrix(adjA,n,n)
Y2M=matrix(X2,n,n)
Risk_n=( (AM<=t(Y2M)&t(Y2M)<=Y2M) )

BBW=matrix(NA,n,1+p)
BBW.sd=list()
PIN=list()
double_check_v=rep(NA,n)
nrisk=matrix(NA,n,1)
nYin=matrix(NA,n,1)
INDIN=list()
XIN=list()

for(tt in 1:n){
if(saving=="TRUE"){gc()}
Yin=1*N1t[(Risk_n[,tt])==1,tt]
nin=length(Yin)
nrisk[tt,1]=nin
nYin[tt,1]=sum(Yin)
if(nin>1) Xin=X[(Risk_n[,tt])==1,]
if(nin==1)Xin=matrix(X[(Risk_n[,tt])==1,],nin,p,byrow=TRUE)

Result_glm=glm(Yin~Xin,family="binomial")

if(is.na(sum(Result_glm$coef))==TRUE) BBW[tt,]=rep(0,1+p) else BBW[tt,]=Result_glm$coef 
#BBW[tt,]=Result_glm$coef 

XIN[[tt]]=cbind(1,Xin)

indin=1:n
INDIN[[tt]]=(1:n)[Risk_n[,tt]]

check_v=vcov(Result_glm)
double_check_v[tt]=max(check_v)
if( is.na(sum(check_v))==TRUE)BBW.sd[[tt]]=matrix(0,(p+1),(p+1)) else BBW.sd[[tt]]=check_v     
           }

double_check_v=ifelse(is.na(double_check_v)==TRUE,10^7,double_check_v)
impu0L=which(double_check_v>10^4 )

################################################################################


for( inp in impu0L){
if(saving=="TRUE"){gc()}
BBW[inp,]=(1:(p+1))*0
BBW.sd[[inp ]]=matrix(0,p+1,p+1)
}


iip=0
check_0=apply(BBW,1,sum)
for(ii in 1:n){
if(saving=="TRUE"){gc()}
if(check_0[ii]==0){ iip=iip+1 } 

if(check_0[ii]!=0& iip>0 & ii!=n ){
BBW_proc=matrix(0,iip,p+1)
T_BBW_proct=t(BBW_proc)
T_BBW_proct[,1:iip]=BBW[ii,]
BBW[(ii-iip):(ii-1),]=t(T_BBW_proct)
BBW.sd[(ii-iip):(ii-1)]=BBW.sd[ii]
iip=0

}

if(iip>0&ii==n){
BBW_proc=matrix(0,(iip+1),p+1)
T_BBW_proct=t(BBW_proc)
T_BBW_proct[,1:(iip+1)]=BBW[(ii-iip-1),]
BBW[(ii-iip):ii,]=t(T_BBW_proct)
BBW.sd[(ii-iip):(ii)]=BBW.sd[(ii-iip-1)]
}
if(saving=="TRUE"){gc()}
}

for(pp in 1:n){
if(saving=="TRUE"){gc()}
expXin=exp(apply(t(XIN[[pp]])*BBW[pp,],2,sum))
pin=rep(0,n)
pin[INDIN[[pp]]]= 1/(1+1/expXin)
PIN[[pp]]=pin
}

#####################################################################################
if(saving=="TRUE"){gc()}



if(center=="FALSE")w1_zb1=apply(t(BBW)*c(1,zb1,rep(0,(p-1))),2,sum )
if(center=="TRUE") w1_zb1=apply(t(BBW)*c(1,zb1,apply(as.matrix(X[,-1]),2,mean)),2,sum)
if(center=="SET")  w1_zb1=apply(t(BBW)*c(1,zb1,setting),2,sum )

pw1_zb1=1/(1+1/exp(w1_zb1))
Lct=which(is.na(pw1_zb1)==1|pw1_zb1=="NaN")
pw1_zb1[Lct]=nYin[Lct]/nrisk[Lct]


pw0_zb1=1-pw1_zb1
pw0_zb1=pw0_zb1



if(sum(EFFECT%in%"IE")==1){
if(center=="FALSE")w1_zb2=apply(t(BBW)*c(1,zb2,rep(0,(p-1))),2,sum )
if(center=="TRUE") w1_zb2=apply(t(BBW)*c(1,zb2,apply(as.matrix(X[,-1]),2,mean)),2,sum)
if(center=="SET")  w1_zb2=apply(t(BBW)*c(1,zb2,setting),2,sum )

pw1_zb2=1/(1+1/exp(w1_zb2))
Lct=which(is.na(pw1_zb2)==1|pw1_zb1=="NaN")
pw1_zb2[Lct]=nYin[Lct]/nrisk[Lct]

pw0_zb2=1-pw1_zb2
pw0_zb2=pw0_zb2
                 }

#BBW=ifelse(is.na(BBW)=="TRUE",0,BBW)
#BBW.sd=ifelse(is.na(BBW.sd)=="TRUE",0,BBW.sd)


if(sum(EFFECT%in%"DE")==1){
dH2z_21=DL0X_za2*pw0_zb1+DL1X_za2*pw1_zb1
dH2z_11=DL0X_za1*pw0_zb1+DL1X_za1*pw1_zb1
DE=cumsum( (dH2z_21-dH2z_11))

}

if(sum(EFFECT%in%"IE")==1){
dH2z_22=DL0X_za2*pw0_zb2+DL1X_za2*pw1_zb2
dH2z_21=DL0X_za2*pw0_zb1+DL1X_za2*pw1_zb1
IE=cumsum( (dH2z_22-dH2z_21) )

}



#########################################################################################
#Mij
#########################################################################################
if(saving=="TRUE"){gc()}
all.alpha=list()
Mij=matrix(0,(p+1)*n,(p+1)*n)

for(ij in 1 :(n-1)){
if(saving=="TRUE"){gc()}
alphij=function(x){
cand=intersect(INDIN[[ij]],INDIN[[x]])
Win=PIN[[ij]]*(1-PIN[[ij]])*(1-PIN[[x]])
Win=Win[cand]
LL=length(Win[Win!=0])
   if(D[x]>0&D[ij]>0){
if(LL>1)  W=diag(Win) else W=as.matrix(Win)
if(LL>1)  Xinproc=cbind(1,X[cand,]) else Xinproc=t(as.matrix(c(1,X[cand,])))
if(LL>1)  ans = BBW.sd[[ij]] %*% t(Xinproc) %*% W %*% Xinproc%*% BBW.sd[[x]] 
if(LL==1) ans = BBW.sd[[ij]] %*% (t(Xinproc) * Win)  %*% Xinproc%*% BBW.sd[[x]]
if(LL==0) ans=matrix(0,p+1,p+1)
}
if(D[x]==0|D[ij]==0) ans=matrix(0,p+1,p+1)
 ans                    }
#--------------------------------------------
 if(D[ij]==1){
List.alpha=apply(as.matrix((ij+1):n),1,alphij)
subMij=matrix(List.alpha,(p+1)*(n-ij),p+1,byrow=TRUE)
L_subMij=NROW(subMij)
Mij[ (n*(p+1)-L_subMij+1) : ((p+1)*n ), ( (ij-1)*(p+1)+1) : (ij*(p+1) )]=subMij
 }
}

Mij=Mij+t(Mij)

for(reij in 1:n){
if(saving=="TRUE"){gc()}
   if(D[reij ]==1){
Mij[ (1:(p+1))+(reij-1)*(p+1) ,(1:(p+1))+(reij-1)*(p+1) ]=BBW.sd[[reij ]]
   }
}


#Det_M=which(apply(Mij,1,sum)==0)
#if(length(Det_M)>0){
#Mij=Mij[-Det_M,-Det_M]}
##########################################################################################
# MB1
##########################################################################################
if(saving=="TRUE"){gc()}
J1.11=J1.12=J1.21=J1.22=J1.13=J1.31=J1.33=J1.32=J1.23=NULL

pdR_F=function(x){
LOC=sum(x>=c(Y21)) 
if(LOC>0) sum(diff_dR1[1:LOC]) else 0
}
pdR_T1=apply(matrix(AY11),1,pdR_F)
EXPB=exp(apply(t(X21)*B21$coef,2,sum))
diffY21_Y11=(cumsum(diff_dR1)-pdR_T1)

J1.11=sum( X21[,1]^2*EXPB* diffY21_Y11  )
##
diff_Risk_n1_jump=diff_Risk_n1[,diff_dR1>0]
if( length(as.matrix(X[,-1]))!=0) {
J1.12=J1.21=apply(as.matrix(X21[,-1]) * (X21[,1]*EXPB* diffY21_Y11 ) , 2 , sum)
J1.12=t(J1.12) 
##
J1.22=t( as.matrix(X21[,-1]) * (EXPB * diffY21_Y11 )^0.5 )%*%( as.matrix(X21[,-1]) * ( EXPB* diffY21_Y11 )^0.5 )
##
J1.23=NULL
for(j23 in 1: NCOL(diff_Risk_n1_jump) ){
J1.23=cbind(J1.23, apply(as.matrix(X21[,-1])* (diff_Risk_n1_jump[,j23]*EXPB),2,sum))
}
J1.32=J1.23
J1.32=t(J1.32)
}
##
J1.13=J1.31=apply( diff_Risk_n1_jump * (X21[,1] * EXPB) , 2 , sum )
J1.31=as.matrix(J1.31)
##
J1.33=( diag(1/diff_dR1[diff_dR1>0]) )^2
##
MB1=cbind(
      c(J1.11,J1.21,J1.31),
      rbind(J1.12,J1.22,J1.32),
      rbind(J1.13,J1.23,J1.33)
)
#MB1=solve(MB1/sum(dN21))
MB1=solve(MB1)


##########################################################################################
# MB0
##########################################################################################
if(saving=="TRUE"){gc()}
J0.11=J0.12=J0.21=J0.22=J0.13=J0.31=J0.33=J0.32=J0.23=NULL
MB0=NULL
if(sum(B20$coef%in%NA)==0&length(B20)!=0&sum(workL0)>0){

pdR_F_0=function(x){
LOC=sum(x>=c(Y20)) 
if(LOC>0) sum(diff_dR0[1:LOC]) else 0
}
pdR_T1_0=apply(matrix(Y10),1,pdR_F_0)*0
EXPB_0=exp(apply(t(X20)*B20$coef,2,sum))
diffY20_Y10=(cumsum(diff_dR0)-pdR_T1_0)
diffY20_Y10=ifelse(is.na(diffY20_Y10)=="TRUE",0,diffY20_Y10)

J0.11=sum(X20[,1]^2*EXPB_0* diffY20_Y10)
##
diff_Risk_n0_jump=diff_Risk_n0[,diff_dR0>0]
if( length(as.matrix(X[,-1]))!=0) {
J0.12=J0.21=apply(as.matrix(X20[,-1])* (X20[,1]*EXPB_0*diffY20_Y10),2,sum)
J0.12=t(J0.12)
##
J0.22=t(as.matrix(X20[,-1]) * (EXPB_0*diffY20_Y10 )^0.5 )%*%( as.matrix(X20[,-1]) * ( EXPB_0* diffY20_Y10 )^0.5 )
##
J0.23=NULL
for(j23 in 1: NCOL(diff_Risk_n0_jump) ){
J0.23=cbind(J0.23, apply(as.matrix(X20[,-1])* (diff_Risk_n0_jump[,j23]*EXPB_0),2,sum))
}
J0.32=J0.23
J0.32=t(J0.32)
}
##

J0.13=J0.31=apply( diff_Risk_n0_jump * (X20[,1]*EXPB_0),2,sum)
J0.31=as.matrix(J0.31)
##
J0.33=( diag(1/diff_dR0[diff_dR0>0]) )^2
MB0=cbind(
      c(J0.11,J0.21,J0.31),
      rbind(J0.12,J0.22,J0.32),
      rbind(J0.13,J0.23,J0.33)
)
#MB0=solve(MB0/sum(1-dN21))
MB0=solve(MB0)
}
##########################################################################################
if(saving=="TRUE"){gc()}
imput_for_a_1=function(x){
TOT_1=sum(Y21[D21>0]<=x)
if(TOT_1==0) ans1=0
if(TOT_1>0&TOT_1<=sum(dN21==1)) ans1= TOT_1
if(TOT_1>sum(dN21==1)) ans1= sum(dN21==1)
ans1
}


imput_for_a_0=function(x){
TOT_0=sum(Y20[D20>0]<=x)
if(TOT_0==0) ans0=0
if(TOT_0>0&TOT_0<=sum(workL0)) ans0= TOT_0
if(TOT_0>sum(workL0)) ans0= sum(workL0)
ans0
}

imput_for_X2=function(x){
TOT_X2=sum(X2[D>0]<=x)
if(TOT_X2==0) ansX2=0
if(TOT_X2>0&TOT_X2<=sum(D==1)) ansX2= TOT_X2
if(TOT_X2>sum(D==1)) ansX2= sum(D==1)
ansX2
}


impaz1=apply(matrix(X2),1,imput_for_a_1)
impaz0=apply(matrix(X2),1,imput_for_a_0)
#impX2 =apply(matrix(X2),1,imput_for_X2)


check_in_F=function(imp_v,vor){
#imp_v=impro;vor=mmpro

if( NCOL(vor)==1) ans=c(rep(0,sum(imp_v==0) ),vor[imp_v] ) else ans=rbind( matrix(0,sum(imp_v==0),NCOL(vor) ),vor[imp_v,] )
ans
}
##########################################################################################
#nn1=n/sum(dN21==1)
#Nn1=nn1^0.5
#nn0=Nn0=sum(workL0)
#nn0=(n/ifelse(Nn0==0,1,Nn0))^0.5
#Nn0=nn0^0.5
Nn1=Nn0=nn1=nn0=1
if(saving=="TRUE"){gc()}
##########################################################################################
#part.L.az
part.L.az=function(za,zb,n1){
#za=za2;zb=zb2;center="TRUE";n1=0

if(n1==1)BB=B21$coef else BB=B20$coef
if(n1==1)D=D21 else D=D20
if(n1==1)diff_dR=diff_dR1 else diff_dR=diff_dR0



if(center=="FALSE"){ XBW=apply(t(BBW)*c(1,zb,rep(0,(p-1))),2,sum )                   ;XBL=sum(BB*c(za,rep(0,(p-1))) )                  }
if(center=="TRUE") { XBW=apply(t(BBW)*c(1,zb,apply(as.matrix(X[,-1]),2,mean)),2,sum) ;XBL=sum(BB*c(za,apply(as.matrix(X[,-1]),2,mean)))}
if(center=="SET")  { XBW=apply(t(BBW)*c(1,zb,setting),2,sum )                        ;XBL=sum(BB*c(za,setting))                        }

if(n1==1)XBW=XBW[dN21==n1]
if(n1==0)XBW=XBW[workL0]

ans=matrix(0,n,1)

org=D*diff_dR* (-1)^(n1+1)*zb*exp(XBL+XBW) /( (1+exp(XBW) )^2)

if(n1==1)ans[dN21==n1,]=org
if(n1==0)ans[workL0,]=org

ifelse(is.na(ans)==TRUE,0,ans)
}
##########################################################################################
#part.L.ax
if(saving=="TRUE"){gc()}
part.L.ax=function(za,zb,n1){
#za=za2;zb=zb2;center="TRUE";n1=0
if(n1==1)BB=B21$coef else BB=B20$coef
if(n1==1)D=D21 else D=D20
if(n1==1)diff_dR=diff_dR1 else diff_dR=diff_dR0


if(center=="FALSE"){ XBW=apply(t(BBW)*c(1,zb,rep(0,(p-1))),2,sum )                   ;XBL=sum(BB*c(za,rep(0,(p-1))) )                  ; X_part=c(1,0) }
if(center=="TRUE") { XBW=apply(t(BBW)*c(1,zb,apply(as.matrix(X[,-1]),2,mean)),2,sum) ;XBL=sum(BB*c(za,apply(as.matrix(X[,-1]),2,mean))); X_part=c(1,apply(as.matrix(X[,-1]),2,mean)) }
if(center=="SET")  { XBW=apply(t(BBW)*c(1,zb,setting),2,sum )                        ;XBL=sum(BB*c(za,setting))                        ; X_part=c(1,setting) }

if(n1==1)XBW=XBW[dN21==n1]
if(n1==0)XBW=XBW[workL0]

org=(D*diff_dR* (-1)^(n1+1)*exp(XBL+XBW) /( (1+exp(XBW) )^2))%x% t(X_part)

ans=matrix(0,n,p)


if(n1==1)ans[dN21==n1,]=org
if(n1==0)ans[workL0,]=org

ifelse(is.na(ans)==TRUE,0,ans)


}

##########################################################################################
#part.L.bz
if(saving=="TRUE"){gc()}
part.L.bz=function(za,zb,n1){
#za=za1;zb=zb2;center="TRUE";n1=0
if(n1==1)BB=B21$coef else BB=B20$coef
if(n1==1)D=D21 else D=D20
if(n1==1)diff_dR=diff_dR1 else diff_dR=diff_dR0

if(center=="FALSE"){ XBW=apply(t(BBW)*c(1,zb,rep(0,(p-1))),2,sum )                   ;XBL=sum(BB*c(za,rep(0,(p-1))) )                  ; X_part=c(1,0) }
if(center=="TRUE") { XBW=apply(t(BBW)*c(1,zb,apply(as.matrix(X[,-1]),2,mean)),2,sum) ;XBL=sum(BB*c(za,apply(as.matrix(X[,-1]),2,mean))); X_part=c(1,apply(as.matrix(X[,-1]),2,mean)) }
if(center=="SET")  { XBW=apply(t(BBW)*c(1,zb,setting),2,sum )                        ;XBL=sum(BB*c(za,setting))                        ; X_part=c(1,setting) }
if(n1==1)XBW=XBW[dN21==n1]
if(n1==0)XBW=XBW[workL0]

ans=za*exp(XBL)* cumsum( (exp(XBW)^n1) / ( 1+exp(XBW) ) *diff_dR)[D>0]

ifelse(is.na(ans)==TRUE,0,ans)
}
##########################################################################################
#part.L.bx
if(saving=="TRUE"){gc()}
part.L.bx=function(za,zb,n1){
#za=za1;zb=zb2;center="TRUE";n1=0
if(n1==1)BB=B21$coef else BB=B20$coef
if(n1==1)D=D21 else D=D20
if(n1==1)diff_dR=diff_dR1 else diff_dR=diff_dR0

if(center=="FALSE"){ XBW=apply(t(BBW)*c(1,zb,rep(0,(p-1))),2,sum )                   ;XBL=sum(BB*c(za,rep(0,(p-1))) )                  ; X_part=c(1,0) }
if(center=="TRUE") { XBW=apply(t(BBW)*c(1,zb,apply(as.matrix(X[,-1]),2,mean)),2,sum) ;XBL=sum(BB*c(za,apply(as.matrix(X[,-1]),2,mean))); X_part=c(1,apply(as.matrix(X[,-1]),2,mean)) }
if(center=="SET")  { XBW=apply(t(BBW)*c(1,zb,setting),2,sum )                        ;XBL=sum(BB*c(za,setting))                        ; X_part=c(1,setting) }

if(n1==1)XBW=XBW[dN21==n1]
if(n1==0)XBW=XBW[workL0]


org=(cumsum(exp(XBL)*(exp(XBW)^n1) / ( 1+exp(XBW) ) *diff_dR) %x% t(X_part[-1]))
ans=as.matrix(org[D>0,])

ifelse(is.na(ans)==TRUE,0,ans)
}



##########################################################################################
#part.L.L
if(saving=="TRUE"){gc()}
part.L.L=function(za,zb,n1){
#za=za1;zb=zb2;center="TRUE";n1=1
if(n1==1)BB=B21$coef else BB=B20$coef
if(n1==1)D=D21 else D=D20

if(center=="FALSE"){ XBW=apply(t(BBW)*c(1,zb,rep(0,(p-1))),2,sum )                   ;XBL=sum(BB*c(za,rep(0,(p-1))) )                  ; X_part=c(1,0) }
if(center=="TRUE") { XBW=apply(t(BBW)*c(1,zb,apply(as.matrix(X[,-1]),2,mean)),2,sum) ;XBL=sum(BB*c(za,apply(as.matrix(X[,-1]),2,mean))); X_part=c(1,apply(as.matrix(X[,-1]),2,mean)) }
if(center=="SET")  { XBW=apply(t(BBW)*c(1,zb,setting),2,sum )                        ;XBL=sum(BB*c(za,setting))                        ; X_part=c(1,setting) }

if(n1==1)XBW=XBW[dN21==n1]
if(n1==0)XBW=XBW[workL0]



ans=(exp(XBL)*D* (exp(XBW)^n1) / ( 1+exp(XBW) ))[D>0] 

ifelse(is.na(ans)==TRUE,0,ans)
}

#######################################################################################
#DE-part 
######################################################################################
if(saving=="TRUE"){gc()}
if(sum(EFFECT%in%c("DE"))==1){
#D_az=check_in_F(impaz1,as.matrix(part.L.az(za2,zb1,1)/Nn1-part.L.az(za1,zb1,1)/Nn1))+ check_in_F(impaz0,as.matrix(part.L.az(za2,zb1,0)/Nn0-part.L.az(za1,zb1,0)/Nn0))
#D_az=as.matrix(D_az)
#D_ax=check_in_F(impaz1,as.matrix(part.L.ax(za2,zb1,1)/Nn1-part.L.ax(za1,zb1,1)/Nn1))+ check_in_F(impaz0,as.matrix(part.L.ax(za2,zb1,0)/Nn0-part.L.ax(za1,zb1,0)/Nn0))
#D_ax=as.matrix(D_ax)

D_az=part.L.az(za2,zb1,1)/Nn1-part.L.az(za1,zb1,1)/Nn1+ part.L.az(za2,zb1,0)/Nn0-part.L.az(za1,zb1,0)/Nn0
D_az=as.matrix(D_az)
D_ax=part.L.ax(za2,zb1,1)/Nn1-part.L.ax(za1,zb1,1)/Nn1+ part.L.ax(za2,zb1,0)/Nn0-part.L.ax(za1,zb1,0)/Nn0
D_ax=as.matrix(D_ax)

D_bz=as.matrix((part.L.bz(za2,zb1,1)-part.L.bz(za1,zb1,1)))


D_bx=NULL
if( length(as.matrix(X[,-1]))!=0) {
D_bx=as.matrix((part.L.bx(za2,zb1,1)-part.L.bx(za1,zb1,1)))}
D_bL=as.matrix((part.L.L(za2,zb1,1)-part.L.L(za1,zb1,1)))


D_bz_0=as.matrix((part.L.bz(za2,zb1,0)-part.L.bz(za1,zb1,0)))
D_bx_0=NULL
if( length(as.matrix(X[,-1]))!=0) {
D_bx_0=as.matrix((part.L.bx(za2,zb1,0)-part.L.bx(za1,zb1,0)))
}
D_bL_0=as.matrix((part.L.L(za2,zb1,0)-part.L.L(za1,zb1,0)))

DE_alpha_sd=rep(NA,n)
for(aa in 1:n){
if(saving=="TRUE"){gc()}
if(aa==1){DE_all_alpha=c(D_ax[1:aa,1],D_az[1:aa] , D_ax[1:aa,-1])} else { DE_all_alpha= c(t(cbind(D_ax[1:aa,1],D_az[1:aa] , D_ax[1:aa,-1]))) }
DE_alpha_sd[aa]=t(DE_all_alpha)%*%Mij[1:(aa*(p+1)),1:(aa*(p+1))]%*%(DE_all_alpha)
}


D_Ln1=which(dN21==1&D>0)
D_n_ln1=length(D_Ln1)

D_Ln0=which(D20[workL0]>0)
D_n_ln0=length(D_Ln0)


DE_cox_sd=rep(0,D_n_ln1)
DE_cox_sd_0=rep(0,D_n_ln0)
reDE_cox_jump=reDE_cox_jump_0=0


for(cc in 1:D_n_ln1){
if(saving=="TRUE"){gc()}
DE_cox_sd[cc ]=t(c(D_bz[cc],D_bx[cc,],D_bL[1:cc,]))%*%MB1[1:(p+cc),1:(p+cc) ]%*%c(D_bz[cc],D_bx[cc,],D_bL[1:cc,])
}

reDE_cox_jump=check_in_F(impaz1,as.matrix(DE_cox_sd))/nn1



if(D_n_ln0>0){
if(saving=="TRUE"){gc()}
for(cc in 1:D_n_ln0){
DE_cox_sd_0[cc]=t(c(D_bz_0[cc],D_bx_0[cc,],D_bL_0[1:cc,]))%*%MB0[1:(p+cc),1:(p+cc) ]%*%c(D_bz_0[cc],D_bx_0[cc,],D_bL_0[1:cc,])
}
reDE_cox_jump_0=check_in_F(impaz0,as.matrix(DE_cox_sd_0))/nn0
}



DE_cox_sd=(reDE_cox_jump+reDE_cox_jump_0)
DE.sd=(DE_cox_sd+DE_alpha_sd)
}
#############################################################################################
#IE-part 
#############################################################################################
if(sum(EFFECT%in%c("IE"))==1){

#I_az=check_in_F(impaz1,as.matrix(part.L.az(za2,zb2,1)/Nn1-part.L.az(za2,zb1,1)/Nn1))+ check_in_F(impaz0,as.matrix(part.L.az(za2,zb2,0)/Nn0-part.L.az(za2,zb1,0)/Nn0))
#I_az=as.matrix(I_az)
#I_ax=check_in_F(impaz1,as.matrix(part.L.ax(za2,zb2,1)/Nn1-part.L.ax(za2,zb1,1)/Nn1))+ check_in_F(impaz0,as.matrix(part.L.ax(za2,zb2,0)/Nn0-part.L.ax(za2,zb1,0)/Nn0))
#I_ax=as.matrix(I_ax)
#I_bz=as.matrix((part.L.bz(za2,zb2,1)-part.L.bz(za2,zb1,1)))

I_az=part.L.az(za2,zb2,1)/Nn1-part.L.az(za2,zb1,1)/Nn1+ part.L.az(za2,zb2,0)/Nn0-part.L.az(za2,zb1,0)/Nn0
I_az=as.matrix(I_az)
I_ax=part.L.ax(za2,zb2,1)/Nn1-part.L.ax(za2,zb1,1)/Nn1+ part.L.ax(za2,zb2,0)/Nn0-part.L.ax(za2,zb1,0)/Nn0
I_ax=as.matrix(I_ax)

I_bz=as.matrix((part.L.bz(za2,zb2,1)-part.L.bz(za2,zb1,1)))


I_bx=NULL
if( length(as.matrix(X[,-1]))!=0) {
I_bx=as.matrix((part.L.bx(za2,zb2,1)-part.L.bx(za2,zb1,1)))}

I_bL=as.matrix((part.L.L(za2,zb2,1)-part.L.L(za2,zb1,1)))

I_bz_0=as.matrix((part.L.bz(za2,zb2,0)-part.L.bz(za2,zb1,0)))

I_bx_0=NULL
if( length(as.matrix(X[,-1]))!=0) {
I_bx_0=as.matrix((part.L.bx(za2,zb2,0)-part.L.bx(za2,zb1,0)))
}
I_bL_0=as.matrix((part.L.L(za2,zb2,0)-part.L.L(za2,zb1,0)))


IE_alpha_sd=rep(NA,n)
for(aa in 1:n){
if(saving=="TRUE"){gc()}
if(aa==1){IE_all_alpha=c(        I_ax[1:aa,1],I_az[1:aa], I_ax[1:aa,-1])} else {
         IE_all_alpha= c(t(cbind(I_ax[1:aa,1],I_az[1:aa], I_ax[1:aa,-1]))) }
IE_alpha_sd[aa]=t(IE_all_alpha)%*%Mij[1:(aa*(p+1)),1:(aa*(p+1))]%*%(IE_all_alpha)
}



I_Ln1=which(dN21==1&D>0)
I_n_ln1=length(I_Ln1)
DD=D
DD[dN21==1]=0
I_Ln0=which(DD[workL0]>0)
I_n_ln0=length(I_Ln0)


IE_cox_sd=rep(0,I_n_ln1)
IE_cox_sd_0=rep(0,I_n_ln0)
reIE_cox_jump=reIE_cox_jump_0=0


for(cc in 1:I_n_ln1){
if(saving=="TRUE"){gc()}
IE_cox_sd[cc ]=t(c(I_bz[cc],I_bx[cc,],I_bL[1:cc,]))%*%MB1[1:(p+cc),1:(p+cc) ]%*%c(I_bz[cc],I_bx[cc,],I_bL[1:cc,])
}
reIE_cox_jump=check_in_F(impaz1,as.matrix(IE_cox_sd))/nn1



if(I_n_ln0>0){
for(cc in 1:I_n_ln0){
if(saving=="TRUE"){gc()}
IE_cox_sd_0[cc]=t(c(I_bz_0[cc],I_bx_0[cc,],I_bL_0[1:cc,]))%*%MB0[1:(p+cc),1:(p+cc) ]%*%c(I_bz_0[cc],I_bx_0[cc,],I_bL_0[1:cc,])
}
reIE_cox_jump_0=check_in_F(impaz0,as.matrix(IE_cox_sd_0))/nn0
}



IE_cox_sd=(reIE_cox_jump+reIE_cox_jump_0)
IE.sd=(IE_cox_sd+IE_alpha_sd)

}

############################################################################################################
#time=t
#cbind(BBW.sd,BBW,pw1_zb1,pw1_zb2)

LT_F=function(x){
LT=sum(x>X2)
if(length(LT)==0) ans=0
if(LT>=0) ans=LT
if(LT>n) ans=n
ans
}

LT=apply(matrix(time),1,LT_F)
T_N0=sum(LT==0)
L_im=LT[LT>0]
if(sum(EFFECT%in%c("DE"))==1){
DE.out=c(rep(0,T_N0),DE[L_im])
DE.sd.out=c(rep(0,T_N0),DE.sd[L_im])
}
if(saving=="TRUE"){gc()}
if(sum(EFFECT%in%c("IE"))==1){
IE.out=c(rep(0,T_N0),IE[L_im])
IE.sd.out=c(rep(0,T_N0),IE.sd[L_im])
}
#####################################

P.F=function(ptime,obstime,dFt){
#ptime=time; obstime=X2; dFt=H2z0_za2;t=ptime[11]
pred.f=function(t){
Ft=cumsum(dFt)
Cand=sum(t>=obstime)
if(Cand==0)ans=0
if(Cand>0)ans=Ft[Cand]
ans
}

apply(matrix(ptime),1,pred.f)

}
#####################################
LA_Sz0_za2=NULL
LA_Sz1_za2=NULL
LA_Sz0_za1=NULL
LA_Sz1_za1=NULL

#####################################
if(sum(EFFECT%in%c("IE","DE"))==2){
lambda_Sz0_za2=DL0X_za2
lambda_Sz1_za2=DL1X_za2
lambda_Sz0_za1=DL0X_za1
lambda_Sz1_za1=DL1X_za1

LA_Sz0_za2=P.F(time,X2,DL0X_za2)
LA_Sz1_za2=P.F(time,X2,DL1X_za2)
LA_Sz0_za1=P.F(time,X2,DL0X_za1)
LA_Sz1_za1=P.F(time,X2,DL1X_za1)
}
#####################################
report=list(
DE=DE.out,
IE=IE.out,
DE.sd=DE.sd.out,
IE.sd=IE.sd.out,
LA_Sz0_za2=LA_Sz0_za2,
LA_Sz1_za2=LA_Sz1_za2,
LA_Sz0_za1=LA_Sz0_za1,
LA_Sz1_za1=LA_Sz1_za1
)
return(report)
}
