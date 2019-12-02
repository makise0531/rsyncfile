#include "TMath.h"
//#include "myinclude.h"

Double_t test(Double_t a){
        cout << TMath::Factorial(a) << endl;
        return 0;
}

Double_t ClebschGordan(Double_t j1in, Double_t m1in, Double_t j2in, Double_t m2in, Double_t jin, Double_t min){
        Double_t cg=0;

        Double_t j1,m1,j2,m2,j,m;
        j1=j1in;   m1=m1in;   j2=j2in;    m2=m2in;    j=jin;    m=min;

        Int_t kmin[3]={};
        Int_t kmax[3]={};

        kmin[0]=0; kmin[1]=j2-j-m1;  kmin[2]=j1+m2-j;
        kmax[0]=j1-m1; kmax[1]=j2+m2;  kmax[2]=j1+j2-j;

        Int_t kmax2=100;
        Int_t kmin2=0.;

        for(Int_t i=0; i<3; i++) {
                //    cout << "   kmin " << kmin[i] << "  kmax " << kmax[i] << endl;
                if(kmax2>kmax[i]) { //　配列内での最小値を求める
                        kmax2=kmax[i];
                }
                if(kmin2<kmin[i]) { //　配列内での最大値を求める
                        kmin2=kmin[i];
                }
        }

        //  cout << kmin2 << "   " << kmax2 << endl;

        Double_t ksum=0;

        if(m1+m2==m) {
                cg = sqrt(2*j+1);
                cg *= sqrt(TMath::Factorial(j1+j2-j));
                cg *= sqrt(TMath::Factorial(j2+j-j1));
                cg *= sqrt(TMath::Factorial(j+j1-j2));
                cg *= sqrt(TMath::Factorial(j1+m1));
                cg *= sqrt(TMath::Factorial(j1-m1));
                cg *= sqrt(TMath::Factorial(j2+m2));
                cg *= sqrt(TMath::Factorial(j2-m2));
                cg *= sqrt(TMath::Factorial(j+m));
                cg *= sqrt(TMath::Factorial(j-m));
                cg /= sqrt(TMath::Factorial(j+j1+j2+1));

                for(Int_t i=kmin2; i<=kmax2; i++) {
                        ksum += pow(-1,i)/(TMath::Factorial(i)*TMath::Factorial(j1-m1-i)*TMath::Factorial(j2+m2-i)*TMath::Factorial(j-j2+m1+i)*TMath::Factorial(j-j1-m2+i)*TMath::Factorial(j1+j2-j-i));
                }
                cg = cg*ksum;
        }

        //  cout << cg << endl;

        return cg;
}


Double_t Racah(Double_t a, Double_t b, Double_t c){

        Double_t rval=0;

        rval = sqrt(TMath::Factorial(a+b-c));
        rval *= sqrt(TMath::Factorial(b+c-a));
        rval *= sqrt(TMath::Factorial(c+a-b));
        rval /= sqrt(TMath::Factorial(a+b+c+1));

        return rval;
}

Double_t Wigner6j(Double_t ein, Double_t ain, Double_t fin, Double_t bin, Double_t din, Double_t cin){
        Double_t Racah(Double_t, Double_t, Double_t);
        Double_t w6j=0;

        Double_t e,a,f,b,d,c;
        e=ein;   a=ain;   f=fin;   b=bin;   d=din;    c=cin;

        Int_t zmax[4]={};
        zmax[0]=2*b;  zmax[1]=b+c-e;  zmax[2]=-a+b+c;  zmax[3]=b-d+f;

        Int_t zmax2=100;

        for(Int_t i=0; i<4; i++) {
                //    cout << zmax[i] << endl;
                if(zmax2>zmax[i]) { //　配列内での最小値を求める
                        zmax2=zmax[i];
                }
        }

        Double_t zsum=0;

        w6j = pow(-1,b+c+e+f);
        w6j *= Racah(a,b,c);
        w6j *= Racah(a,e,f);
        w6j *= Racah(c,d,e);
        w6j *= Racah(b,d,f);
        w6j *= TMath::Factorial(a+b+c+1);
        w6j *= TMath::Factorial(b+d+f+1);
        w6j /= TMath::Factorial(c-d+e);
        w6j /= TMath::Factorial(c+d-e);
        w6j /= TMath::Factorial(a-e+f);
        w6j /= TMath::Factorial(-a+e+f);
        w6j /= TMath::Factorial(b+d-f);

        for(Int_t i=0; i<=zmax2; i++) {
                zsum += pow(-1,i)*TMath::Factorial(2*b-i)*TMath::Factorial(b+c-e+f-i)*TMath::Factorial(b+c+e+f+1-i) / (TMath::Factorial(i)*TMath::Factorial(-a+b+c-i)*TMath::Factorial(b-d+f-i)*TMath::Factorial(a+b+c+1-i)*TMath::Factorial(b+d+f+1-i));
        }

        w6j = w6j*zsum;

        return w6j;

}
/*
   Double_t Wigner9j(Double_t j1in, Double_t j2in, Double_t j3in, Double_t j4in, Double_t j5in, Double_t j6in, Double_t j7in, Double_t j8in, Double_t j9in){

   Double_t Wigner6j(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
   Double_t w9j=0;

   Double_t j1,j2,j3,j4,j5,j6,j7,j8,j9;
   j1=j1in;  j2=j2in;  j3=j3in;  j4=j4in;  j5=j5in;  j6=j6in;  j7=j7in;  j8=j8in;  j9=j9in;

   Int_t xmin[3]={};
   Int_t xmax[3]={};

   xmin[0]=fabs(j1-j4); xmin[1]=fabs(j2-j8);  xmin[2]=fabs(j6-j9);
   xmax[0]=j1+j4; xmax[1]=j2+j8;  xmax[2]=j6+j9;

   Int_t xmax2=100;
   Int_t xmin2=0.;

   for(Int_t i=0;i<3;i++){
    cout << "   xmin " << xmin[i] << "  xmax " << xmax[i] << endl;
    if(xmax2>xmax[i]){  //　配列内での最小値を求める
      xmax2=xmax[i];
    }
    if(xmin2<xmin[i]){  //　配列内での最大値を求める
      xmin2=xmin[i];
    }
   }

   cout << xmin2 << "   " << xmax2 << endl;

   for(Int_t i=xmin2;i<=xmax2;i++){
    w9j = pow(-1,2*i)*(2*i+1);
    cout << w9j << endl;
    w9j *= Wigner6j(j1,j4,j7,j8,j9,i);
    cout << w9j << endl;
    w9j *= Wigner6j(j2,j5,j8,j4,i,j6);
    cout << w9j << endl;
    w9j *= Wigner6j(j3,j6,j9,i,j1,j2);
   }

   cout << w9j << endl;

   return w9j;
   }
 */

Double_t RevTri(Double_t a, Double_t b, Double_t c){

        Double_t rval=0;

        rval = sqrt(TMath::Factorial(a-b+c));
        rval *= sqrt(TMath::Factorial(a+b-c));
        rval *= sqrt(TMath::Factorial(a+b+c+1));
        rval /= sqrt(TMath::Factorial(b+c-a));

        return rval;
}

Double_t Wigner9j(Double_t j1in, Double_t j2in, Double_t j3in, Double_t j4in, Double_t j5in, Double_t j6in, Double_t j7in, Double_t j8in, Double_t j9in){

        Double_t RevTri(Double_t, Double_t, Double_t);
        //Double_t Wigner6j(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t);
        Double_t w9j=0;

        Double_t j1,j2,j3,j4,j5,j6,j7,j8,j9;
        j1=j1in;  j2=j2in;  j3=j3in;  j4=j4in;  j5=j5in;  j6=j6in;  j7=j7in;  j8=j8in;  j9=j9in;

        //  x range
        Int_t xmin[3]={};
        Int_t xmax[3]={};

        xmin[0]=0; xmin[1]=0;  xmin[2]=0;
        xmax[0]=2*j9; xmax[1]=j5-j2+j8;  xmax[2]=j3+j6-j9;

        Int_t xmax2=100;
        Int_t xmin2=0.;

        for(Int_t i=0; i<3; i++) {
                //    cout << "   xmin " << xmin[i] << "  xmax " << xmax[i] << endl;
                if(xmax2>xmax[i]) { //　配列内での最小値を求める
                        xmax2=xmax[i];
                }
                if(xmin2<xmin[i]) { //　配列内での最大値を求める
                        xmin2=xmin[i];
                }
        }
        //  cout << xmin2 << "   " << xmax2 << endl;

        // y range
        Int_t ymin[3]={};
        Int_t ymax[3]={};

        ymin[0]=0;
        ymax[0]=j3-j6+j9;

        Int_t ymax2=100;
        Int_t ymin2=0.;

        for(Int_t i=0; i<1; i++) {
                //    cout << "   ymin " << ymin[i] << "  ymax " << ymax[i] << endl;
                if(ymax2>ymax[i]) { //　配列内での最小値を求める
                        ymax2=ymax[i];
                }
                if(ymin2<ymin[i]) { //　配列内での最大値を求める
                        ymin2=ymin[i];
                }
        }
        //  cout << ymin2 << "   " << ymax2 << endl;

        // z range
        Int_t zmin[2]={};
        //  Int_t zmax[2]={};

        zmin[0]=0; zmin[1]=j1-j4-j7;

        Int_t zmax2=j1-j4+j7;
        Int_t zmin2=0.;

        for(Int_t i=0; i<2; i++) {
                //    cout << "   zmin " << zmin[i] << "  zmax " << zmax[i] << endl;
                //    cout << "   zmin " << zmin[i] << "  zmax " << zmax2 << endl;
                if(zmin2<zmin[i]) { //　配列内での最大値を求める
                        zmin2=zmin[i];
                }
        }
        //  cout << zmin2 << "   " << zmax2 << endl;

        w9j = pow(-1,j7+j8+j9);
        w9j *= RevTri(j2,j1,j3);
        w9j *= RevTri(j4,j5,j6);
        w9j *= RevTri(j9,j3,j6);
        w9j /= RevTri(j2,j5,j8);
        w9j /= RevTri(j4,j1,j7);
        w9j /= RevTri(j9,j7,j8);

        //  cout << zmin2 << "   " << zmax2 << endl;

        Double_t sum=0;
        Double_t sum_temp=0;

        //  cout << zmin2 << "   " << zmax2 << endl;

        for(Int_t i=xmin2; i<=xmax2; i++) { // x
                for(Int_t j=ymin2; j<=ymax2; j++) { // y
                        for(Int_t k=zmin2; k<=zmax2; k++) { // z
                                if(i+k>=j1-j4-j9+j8) {
                                        sum_temp=0;
                                        //	  cout << "aa  " << i << "  " << j << "   " << k << endl;
                                        sum_temp = pow(-1,i+j+k);
                                        sum_temp *= TMath::Factorial(2*j8-i);
                                        sum_temp *= TMath::Factorial(j2+j5-j8+i);
                                        sum_temp /= TMath::Factorial(j5-j2+j8-i);
                                        sum_temp /= TMath::Factorial(j7+j8-j9-i);
                                        sum_temp /= TMath::Factorial(j2-j4+j6-j8+i+j);

                                        sum_temp *= TMath::Factorial(j7-j8+j9+i);
                                        sum_temp *= TMath::Factorial(j5-j4+j6+j);
                                        sum_temp /= TMath::Factorial(j4-j1-j8+j9+i+k);
                                        sum_temp /= TMath::Factorial(j);
                                        sum_temp /= TMath::Factorial(j4+j5-j6-j);

                                        sum_temp *= TMath::Factorial(j3+j6-j9+j);
                                        sum_temp *= TMath::Factorial(j1+j2-j6+j9-j-k);
                                        sum_temp /= TMath::Factorial(j3-j6+j9-j);
                                        sum_temp /= TMath::Factorial(2*j6+1+j);
                                        sum_temp /= TMath::Factorial(k);
                                        sum_temp /= TMath::Factorial(j1+j2-j3-k);

                                        sum_temp *= TMath::Factorial(2*j1-k);
                                        sum_temp *= TMath::Factorial(j4-j1+j7+k);
                                        sum_temp /= TMath::Factorial(j1-j4+j7-k);
                                        sum_temp /= TMath::Factorial(j1+j2+j3+1-k);

                                        sum = sum + sum_temp;
                                        //	  cout << sum << endl;
                                }
                        }
                }
        }

        w9j=w9j*sum;

        /*
           // cal w9j
           for(Int_t i=xmin2;i<=xmax2;i++){
           w9j = pow(-1,2*i)*(2*i+1);
           cout << w9j << endl;
           w9j *= Wigner6j(j1,j4,j7,j8,j9,i);
           cout << w9j << endl;
           w9j *= Wigner6j(j2,j5,j8,j4,i,j6);
           cout << w9j << endl;
           w9j *= Wigner6j(j3,j6,j9,i,j1,j2);
           }
         */
        //  cout << w9j << endl;

        return w9j;
}
