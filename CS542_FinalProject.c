#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    int maxTime = atoi(argv[1]);
    int trackSize = atoi(argv[2]);
    int moveBufCnt = atoi(argv[3]);
    int moveBufSize = (4*moveBufCnt);
    int reqBufSize = (1+3*moveBufCnt);
    
//    //make file
//    FILE * fPtr;
//    char fPath[100];
//    sprintf(fPath,"Problem3/ScatterAllGather/nprocs_%d/ScatterAllGather_N_%d.txt",num_procs,N);
//    int numDoubles = 1 << N;
//    int blockSize = (int) (numDoubles / num_procs);
    
    //walker idx is 2*walker, organized right to left (new walkers are at end of list)
    int* walkers = (int*)calloc(2*trackSize,sizeof(int)); //distToLeft, distToSub
    double* moveBuf = (double*)calloc(moveBufSize,sizeof(double)); //time, walker, move, cleave
    double* outgoingReqBuf = (double*)calloc(reqBufSize,sizeof(double));
    double* incomingReqBuf = (double*)calloc(reqBufSize,sizeof(double));
    double* sendWalker = (double*)calloc(2,sizeof(double)); // send walker buffer
    
    int numWalkers = 0, numWalkersOnSubstrate = 0, numWalkersOnProduct = 0;
    if (rank==0)
    {
        numWalkers = 1;
        numWalkersOnSubstrate=1;
    }
    int leftmost = 0;
    int nu = 0;
    double sigma = 0.0;
    double moveRand = 0.0;// for moving walker
    double time = 0;
    int walker = 0; //moving walker
    int onSubstrate = 1; //bool - for moving walker
    int cleave = 1; //bool - for moving walker
    int rightNeighbor = 0; //for moving walker
    int moveCtr = 0; //moves stored in buffer
    int move = 0; //-1 left, 0 self, 1 right
    int distToEnd = trackSize; //for rightmost walker only
    int distToRight = 0;
    int distToLeft = 0;
    int moveSelf = 0;
    int moveLeft = 0;
    int moveRight = 0;
    int moveSelfOrLeft=0;
    int moveSelfOrRight=0;
    int moveSelfOrLeftOrRight=0;
    int regularMove=0;
    int spawn=0;
    int requestLeft=0;
    int requestRight=0;
    int absorb=0;
    int incomingRequest=0;
    int incomingWalker=0;
    int neighborRank = (rank+1)%1;
    int insertLeft,insertRight;
    int lastMove = 0;
    double requestTime=0.0;
    double incomingWalkerTime=0.0;
    int requestMode=0;
    int send=0;
    int clearMemoryRequest=0;
    int rewindWalker,rewindMove,rewindCleave;
    int accessGranted=0;
    double simulationTimeBacktracked=0;
    int clearMemoryMode=0;
    MPI_Request send_request;
    int trySaveMove=0;
    double clearMemoryTime;
    int clearMemoryGranted;
    int message=0;

    //START TIMING
    MPI_Barrier(MPI_COMM_WORLD);
//    start = MPI_Wtime();
    
    while (time <= maxTime)
    {
        incomingRequest = 0;
        incomingWalker = 0;
        clearMemoryRequest=0;
        clearMemoryGranted=0;
        accessGranted=0;
        
        MPI_Status status;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&message,&status);
        if(message)
            if (status.MPI_TAG==0) incomingRequest=1;
            else if (status.MPI_TAG == 1) incomingWalker=1;
            else if (status.MPI_TAG==3) clearMemoryRequest=1;
        //tag 0 = requestEntry
        //tag 1 = incoming walker
        //tag 2 = accessGranted
        //tag 3 = clearMemoryRequest
        //tag 4 = clearMemoryGranted
        
        if (clearMemoryRequest)
            MPI_Recv(&clearMemoryTime,1,MPI_DOUBLE,neighborRank,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if (time < clearMemoryTime)
            {
                moveCtr=0;
                clearMemoryMode=1;
            }
            else
                MPI_Send(1,1,MPI_INT,neighborRank,4,MPI_COMM_WORLD,send_request); //tag 4 = clearMemoryGranted

        if (incomingRequest)
        {
            MPI_Recv(incomingReqBuf,reqBufSize,MPI_DOUBLE,neighborRank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            lastMove = (int)incomingReqBuf[0];
            requestTime = incomingReqBuf[1+(3*lastMove)];
            if (time < requestTime)
                requestMode=1;
            else if (time > requestTime)
                int i=moveCtr;
                while(outgoingReqBuf[1+3*i]>requestTime)
                    i--;
                if (rank==0)
                    accessGranted=(outgoingReqBuf[1+(3*i)+1]>0);
                else
                    accessGranted=(outgoingReqBuf[1+(3*1)+2]>0);
                MPI_ISend(accessGranted,1,MPI_INT,neighborRank,2,MPI_COMM_WORLD,send_request); //tag 2 = accessGranted
            else
                assert (0);
        }
        
        if (requestMode)
        {
            if (time>requestTime)
            {
                if (rank==0)
                    accessGranted=(distToEnd>0);
                else
                    accessGranted=(walkers[leftmost]>0);
                MPI_ISend(accessGranted,1,MPI_INT,neighborRank,2,MPI_COMM_WORLD,send_request); //tag 2 = accessGranted
                requestMode=0;
            }
        }
        
        insertLeft=0;
        insertRight=0;
        if (incomingWalker)
        {
            MPI_Recv(sendWalker,1,MPI_DOUBLE,neighborRank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            time = sendWalker[0];
            int i=moveCtr;
            while(moveBuf[1+4*i]>time)
            {
                rewindWalker=moveBuf[1+(4*i)+1]
                rewindMove=moveBuf[1+(4*i)+2]
                rewindCleave=moveBuf[1+(4*i)+3]
                walkers[rewindWalker] -= rewindMove;
                walkers[rewindWalker+1] += rewindMove;
                walkers[rewindWalker+1] += rewindCleave;
                i--;
            }
            numWalkers+=1;
            if (rank==0)
                insertRight=1;
            else
                insertLeft=1;
            moveCtr=0;
        }
        
        if (insertRight)
        {
            walkers[0] = walkers[1]-1; // new distToLeft is last rightmost distToSubstrate - 1
            walkers[1] = sendWalker[1]; // new distToSubstrate
            for (i=0;i<(numWalkers-1);i++) // shift all walker indexes to the right (since walker 0 is rightmost)
            {
                walkers[i+2]=walkers[i];
                walkers[i+3]=walkers[i+1];
            }
            distToEnd = 0;
        }
        
        if (insertLeft)
        {
            walkers[numWalkers] = 0;
            walkers[numWalkers+1] = sendWalker[1];
        }
        
        assert (moveCtr <= moveBufCnt);
        double totalRateForSubstrates = r * (double) numWalkersOnSubstrate;
        double totalRateForProducts = (double) numWalkersOnProduct;
        totalRate = totalRateForProducts + totalRateForSubstrates;
    
        double u = uniform01();
        assert (0.0 <= u && u < 1.0);
        rho = totalRate * u;
        
        time += 1.0 / totalRate;

        nu = 0;
        sigma = 0.0;
        while (nu < numWalkers)
        {
            if (walkers [nu] == 0)
                sigma += r;
            else
                sigma += 1.0;
            if (sigma > rho) break;
            nu++;
        }
        
        assert (0 <= nu && nu < numWalkers);
        walker = 2*nu;
        substrate = walker+1;
        leftmost = 2*(numWalkers-1);
        
        //get distances to left and right
        if (walker==0)
            distToRight = distToEnd;
        else
        {
            rightNeighbor = 2*(nu-1);
            distToRight = walkers[rightNeighbor];
        }
        distToLeft = walkers[walker];
        
        //check if on substrate
        onSubstrate=0;
        if (walkers[substrate]==0) onSubstrate=1;
        
        //reset
        moveSelf = 0;
        moveLeft = 0;
        moveRight = 0;
        moveSelfOrLeft=0;
        moveSelfOrRight=0;
        moveSelfOrLeftOrRight=0;
        regularMove=0;
        spawn=0;
        send=0;
        requestLeft=0;
        requestRight=0;
        absorb=0;
            
//----------CLASSIFY WALKER--------------//
        if (walker == leftmost && rank == 0) //origin
        {
            regularMove=1;
            spawn=1;
        }
        else if (walker == leftmost && rank != 0) //leftmost (not origin)
        {
            if (walkers[walker]==0) // walker at interface
            {
                requestLeft=1;
            }
            else
                regularMove=1;
        }
        else if (walker == 0) //rightmost
        {
            if (distToEnd==0)
            {
                if (rank==0)
                {
                    requestRight=1;
                }
                else
                    regularMove=1;
                    absorb=1;
            }
            else
                regularMove=1;
        }
        else //not leftmost or rightmost
            regularMove=1;
        
//-------CHOOSE FEASIBLE POSITION---------//
        
        if (requestLeft)
        {
            if (requestMode) //don't need to request access
            {
                int i=0;
                while (incomingReqBuf[1+(3*i)]<time)
                    i++;
                if (incomingReqBuf[1+3*(i-1)+2]==0) //right gap is zero
                    distToLeft=0;
                else distToLeft=1;
            }
            else
            {
                MPI_Send(outgoingReqBuf,reqBufSize,MPI_DOUBLE,(rank-1),0,MPI_COMM_WORLD);
                MPI_Recv(accessGranted,1,MPI_INT,(rank-1),2,MPI_COMM_WORLD);
                distToLeft=accessGranted;
            }
        }
        
        if (requestRight)
        {
            if (requestMode) //don't need to request access
            {
                int i=0;
                while (incomingReqBuf[1+(3*i)]<time)
                    i++;
                i-=1;
                if (incomingReqBuf[1+(3*i)+1]==0) //left gap is zero
                    distToRight=0;
                else distToRight=1;
            }
            else
            {
                MPI_Send(outgoingReqBuf,reqBufSize,MPI_DOUBLE,(rank+1),0,MPI_COMM_WORLD);
                MPI_Recv(accessGranted,1,MPI_INT,(rank+1),2,MPI_COMM_WORLD);         // WAIT HERE ??********************************************
                distToRight=accessGranted;
            }
        }
        
        if (regularMove)
        {
            if (distToLeft == 0 && distToRight == 0)//blocked on both sides
                moveSelf=1;
            else if (distToLeft>0 && distToRight==0)//blocked only on right
                moveSelfOrLeft=1;
            else if (distToLeft==0 && distToRight>0)//blocked only on left
                moveSelfOrRight=1;
            else if (distToLeft>0 && distToRight>0)//not blocked on either side
                moveSelfOrLeftOrRight=1;
        }
        
        if (moveSelfOrLeft)
        {
            moveRand = uniform01();
            if (moveRand < .5)
                moveSelf=1;
            else
                moveLeft=1;
        }
        
        if (moveSelfOrRight)
        {
            moveRand = uniform01();
            if (moveRand < .5)
                moveSelf=1;
            else
                moveRight=1;
        }
        
        if (moveSelfOrLeftOrRight)
        {
            moveRand = uniform01();
            if (moveRand < 1.0/3.0)
                moveLeft=1;
            else if (moveRand < 2.0/3.0)
                moveSelf=1;
            else
                moveRight=1;
        }
        
//-------------MOVE WALKER---------------//
        
        if (moveSelf)
        {
            move = 0;
            if (onSubstrate==1)
            {
                cleave = 1;
                walkers[substrate] += 1;
                numWalkersOnSubstrate-=1;
                numWalkersOnProduct+=1;
            }
            else cleave = 0;
        }
        
        if (moveRight)
        {
            move = 1;
            walkers[walker] += move; //move right
            if (walker > 0) // right neigbor exists
                walkers[rightNeighbor] -= move;
            if (onSubstrate==0)
            {
                cleave = 0;
                walkers[substrate] -= 1; //distToSubstrate
                if (walkers[substrate] == 0 && !absorb) // stepped onto substrate
                {
                    numWalkersOnSubstrate += 1;
                    numWalkersOnProduct -= 1;
                }
            }
            else
                cleave = 1;
            if (!(absorb || send))
                distToEnd -= 1;
            if (spawn)
            {
                numWalkers += 1;
                numWalkersOnSubstrate += 1;
                walkers[2*numWalkers] = 0;
                walker[2*numWalkers+1] = 0;
            }
            if (requestRight) send=1;
        }
        
        if (moveLeft)
        {
            move = -1;
            walkers[walker] += move; //move left
            if (walker > 0) // right neigbor exists
                walkers[rightNeighbor] -= move;
            else //walker is rightmost
                distToEnd -=move
            if (onSubstrate==0)
            {
                cleave = 0;
                walkers[substrate] -= 1; //distToSubstrate
            }
            else
            {
                cleave = 1;
                walkers[substrate] -= 2*move;
            }
            if (requestLeft) send=1;
        }
        
        if (absorb || (moveRight && send))
        {
            numWalkers-=1;
            distToEnd = distToLeft+1;
            for (i=0;i<(numWalkers-1);i++) // shift all walker indexes to the left (since walker 0 is rightmost)
            {
                walkers[i]=walkers[i+2];
                walkers[i+1]=walkers[i+3];
            }
            if (cleave) numWalkersOnSubstrate -= 1;
            else numWalkersOnProduct -= 1;
        }
        else if (moveLeft && send)
        {
            numWalkers-=1;
            if (cleave) numWalkersOnSubstrate-=1;
            else numWalkersOnProduct -= 1;
        }
        
        if (send)
        {
            sendWalker[0]=time;
            sendWalker[1]=walkers[substrate];
            if requestMode:
            {
                MPI_Isend(0,1,MPI_INT,neighborRank,2,MPI_COMM_WORLD);//send accessDenied since neighbor is waiting for availability
                requestMode=0;
            }
            MPI_Isend(1,1,MPI_INT,neighborRank,4,MPI_COMM_WORLD);//send clearMemoryGranted in case neighbor has requested
            MPI_Isend(sendWalker,1,MPI_DOUBLE,neighborRank,1,MPI_COMM_WORLD,&send_request);
        }

        if (!send) // sending a walker will automatically clear the memory
        {
            if (!clearMemoryMode)
            {
                moveCtr += 1;
                if (moveCtr > moveBufCnt)
                {
                    MPI_Isend(time,1,MPI_DOUBLE,neighborRank,3,MPI_COMM_WORLD,&send_request);
                    MPI_Recv(clearMemoryGranted,1,MPI_INT,neighborRank,4,MPI_COMM_WORLD);
                    moveCtr=0;
                }
                else
                {
                    outgoingReqBuf[0] = (double) moveCtr;
                    outgoingReqBuf[1+3*moveCtr]=time;
                    outgoingReqBuf[1+3*moveCtr+1]=(double) walkers[leftmost];
                    outgoingReqBuf[1+3*moveCtr+2]=(double) distToEnd;
                    //---update history---//
                    moveBuf[1+4*moveCtr] = time;
                    moveBuf[1+4*moveCtr+1] = (double) walker;
                    moveBuf[1+4*moveCtr+2] = (double) move;
                    moveBuf[1+4*moveCtr+3] = (double) cleave;
                }
            }
            else
            {
                if time > clearMemoryTime
                {
                    MPI_Send(1,1,MPI_INT,neighborRank,4,MPI_COMM_WORLD,send_request); //tag 4 = clearMemoryGranted)
                    clearMemoryMode=0;
                }
            }
        }
    }
    MPI_Finalize();
    return (0);
}
 
    
//    if (rank==0)
//    {
//        printf("rank %d sum %e\n",rank,sum);
//        //Save result
//        fPtr = fopen(fPath ,"a");
//        if (fPtr == NULL) exit(EXIT_FAILURE);
//        fprintf(fPtr,"%e\n",time);
//        fclose(fPtr);
//    }
