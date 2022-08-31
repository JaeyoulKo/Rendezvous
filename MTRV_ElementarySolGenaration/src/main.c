#include <stdio.h>
#include "Lambert.h"
#include "parameter.h"

void main(){
    double rm1, rm2, thetadeg, tf;
    double vr1, vt1, vr2, vt2;
    double vr1_c, vt1_c, vr2_c, vt2_c;
    int ind;
    double r1vec[3], v1vec[3], r2vec[3], v2vec[3];
    char buf[5000];
    double tmp1,tmp2,tmp3;
    double Jval[3], sol_gd[10000], x_local[30000];

    printf("hello world!");
    FILE* finp = fopen("lambert_tester.txt", "r");
    for (int i=0 ; i<10 ; i++){
        fgets(buf, 5000, finp);
        fscanf(finp, "%lf,%lf,%lf\t%lf,%lf,%lf\t%lf,%lf,%lf\t%lf,%lf,%lf\t%lf\t%lf\t%lf\t%lf\n" ,
        &r1vec[0],&r1vec[1],&r1vec[2],&v1vec[0],&v1vec[1],&v1vec[2],
        &r2vec[0],&r2vec[1],&r2vec[2],&v2vec[0],&v2vec[1],&v2vec[2],
        &rm1, &rm2, &thetadeg, &tf
        );
        

        FILE* fresult = fopen("lambert_result.txt", "a");
        VLAMB(rm1, rm2, thetadeg , tf, sol_gd);
        fprintf(fresult, "r1vec v1vec r2vec v2vec rm1 rm2 theta(deg) tf \n");
        fprintf(fresult, "%f,%f,%f\t%f,%f,%f\t%f,%f,%f\t%f,%f,%f\t %f \t %f \t %f \t %f \n", 
                r1vec[0],r1vec[1],r1vec[2],v1vec[0],v1vec[1],v1vec[2],
                r2vec[0],r2vec[1],r2vec[2],v2vec[0],v2vec[1],v2vec[2],
                rm1, rm2, thetadeg, tf);
        ind = sol_gd[0];
        for (int i = 0; i < ind; i++){
            vr1 = sol_gd[5 * i + 2];
            fprintf(fresult, "%f\t", sol_gd[5 * i + 2]);
            vt1 = sol_gd[5 * i + 3];
            fprintf(fresult, "%f\t", sol_gd[5 * i + 3]);
            vr2 = sol_gd[5 * i + 4];
            fprintf(fresult, "%f\t", sol_gd[5 * i + 4]);
            vt2 = sol_gd[5 * i + 5];
            fprintf(fresult, "%f\n", sol_gd[5 * i + 5]);
        }
        fclose(fresult);
        fgets(buf, 5000, finp);
        fgets(buf, 5000, finp);
        fscanf(finp, "%lf",&vr1_c,&vt1_c,&vr2_c,&vt2_c);
        fgets(buf,5000,finp);
        printf("%f & %f",vr1_c,vr1);
        if (vr1_c == vr1){
            printf("good!");
        }
        else printf("bad");
        getchar();
    }
}
