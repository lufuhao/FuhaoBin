#include<stdio.h> 
#include<stdlib.h>
size_t maximum=0;
int main(int argc,char *argv[])
{
    void * block;
    void * tmpblock;
    size_t blocksize[]={1024*1024, 1024, 1};
    int i,count;
    for(i=0;i<3;i++){
        for(count=1;;count++){
            block = malloc(maximum+blocksize[i]*count);
            if(block){
                tmpblock = block;
                maximum += blocksize[i]*count;
                free(block);
            }else{
                break;
            }
        }
    }
    printf("maximum malloc size = %lf GB\n",maximum*1.0 / 1024.0 / 1024.0 / 1024.0);
    printf("the address is %x\n",tmpblock);
    printf("the address end is %x\n", tmpblock + maximum);
    //while(1);
}

