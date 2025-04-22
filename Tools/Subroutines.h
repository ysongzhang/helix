#ifndef _Subroutines
#define _Subroutines

/* Defines subroutines for use in both KeyGen and Offline phase
 * Mainly focused around commiting and decommitting to various
 * bits of data
 */

#include "Tools/random.h"
#include "Networking/Player.h"
#include "Tools/Commit.h"

/* Produce a random seed of length len */
void Create_Random_Seed(octet* seed,const PlayerBase& P,int len);



/* Functions to Commit to Seed Values 
 *  This also initialises the PRNG's in G
 */
void Commit_To_Seeds(vector<PRNG>& G,
                     vector< vector<octetStream> >& seeds,
                     vector< vector<octetStream> >& Comm_seeds,
		     vector<octetStream>& Open_seeds,
		     const Player& P,int num_runs);



/* Run just the Commit Protocol for data[i][j] of type T
 *   0 <= i < num_runs
 *   0 <= j < num_players
 * On input data[i][j] is only defined for j=my_number
 */
template<class T>
void Commit(vector< vector<octetStream> >& Comm_data,
            vector<octetStream>& Open_data,
            const vector< vector<T> >& data,const Player& P,int num_runs)
{
  octetStream os;
  int my_number=P.my_num();
  for (int i=0; i<num_runs; i++)
     { os.reset_write_head();
       data[i][my_number].pack(os);
       Comm_data[i].resize(P.num_players());
       Commit(Comm_data[i][my_number],Open_data[i],os,my_number);
       P.Broadcast_Receive(Comm_data[i]);
     }
}


/* Run just the Open Protocol for data[i][j] of type T
 *   0 <= i < num_runs
 *   0 <= j < num_players
 * On output data[i][j] contains all the data
 * If dont!=-1 then dont open this run
 */
template<class T>
void Open(vector< vector<T> >& data,
          const vector< vector<octetStream> >& Comm_data,
          const vector<octetStream>& My_Open_data,
          const Player& P,int num_runs,int dont=-1)
{
  octetStream os;
  int my_number=P.my_num();
  int num_players=P.num_players();
  vector<octetStream> Open_data(num_players);
  for (int i=0; i<num_runs; i++)
     { if (i!=dont)
        { Open_data[my_number]=My_Open_data[i];
          P.Broadcast_Receive(Open_data);
          for (int j=0; j<num_players; j++)
            { if (j!=my_number)
                { if (!Open(os,Comm_data[i][j],Open_data[j],j))
                     { throw invalid_commitment(); }
                  os.reset_read_head();
                  data[i][j].unpack(os);
                }
            }
        }
     }
}


#endif