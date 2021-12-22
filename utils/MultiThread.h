#ifndef __COMMON_MULTITHREAD_MULTITHREAD_H__
#define __COMMON_MULTITHREAD_MULTITHREAD_H__

/**********************************************************
 *
 * @Brief :
 *     A therad pool implementation.
 *
 *     How to use it in main thread :
 *
 *      MultiThread mt;
 *
 *      mt.AddJob(xxx) ....       <-- add all jobs in single thread before all run
 *
 *      mt.Start(8);              <-- 8 worker thread will be created
 *
 *      mt.WaitingStop();         <-- block until all jobs are done.
 *
 *    // Enjob it ~
 * *******************************************************/

#include <functional>
#include <queue>
#include <condition_variable>
#include <mutex>
#include <iostream>
#include <thread>

namespace BGIQD{
    namespace MultiThread{

        // One job --> A function to be called
        typedef std::function<void(void)> Job;

        // A thread safe job pool.
        struct JobQueue
        {
            public:
                JobQueue() {}
                JobQueue( const  JobQueue &) = delete;
                JobQueue & operator=(const JobQueue & ) =delete ;
            public:
                // set in single thread !!!
                void Set( const Job & j)
                {
                    m_queue.push(j);
                }

                std::pair<bool , Job> Get()
                {
                    std::unique_lock<std::mutex> lk(m_mutex);
                    if( m_queue.empty() )
                    {
                        //if( end )
                        return std::make_pair(false , Job());
                    }
                    Job j =  m_queue.front();
                    m_queue.pop();
                    lk.unlock();
                    return std::make_pair(true,j);
                }
            private:
                std::mutex m_mutex;
                std::queue<Job > m_queue;
        };

        // A tireless worker
        struct Thread 
        {
            static void run(JobQueue & queue)
            {
                while(true)
                {
                    auto top = queue.Get();
                    if( top.first == false )
                        break;
                    else
                    {
                        top.second();
                    }
                }
            }
        };

        // The multi-thread manager
        struct MultiThread
        {
            public:

                MultiThread() {}
                MultiThread( const MultiThread & ) = delete ;
                MultiThread& operator=(const MultiThread & ) = delete;

            public:
                // Start multi-thread engine.
                void Start(int thread_num)
                {
                    for( int i = 0 ; i < thread_num ; i++)
                    {
                        m_threads.emplace_back(new std::thread(std::bind(&Thread::run,std::ref(m_queue) )));
                    }
                }

                // After arranging all jobs, the main thread call this to block itself util all jobs are done.
                void WaitingStop()
                {
                    while(true)
                    {
                        size_t ends = 0;
                        for(size_t i = 0 ; i<m_threads.size() ; i++)
                        {
                            if( m_threads[i] )
                            {
                                m_threads[i]->join();
                                delete m_threads[i];
                                m_threads[i] = NULL ; 
                                break;
                            }
                            ends++;
                        }
                        if( ends == m_threads.size() )
                        {
                            break;
                        }
                        else
                        {
                            std::this_thread::sleep_for(std::chrono::milliseconds(1));
                        }
                    }
                }

                // Arrage one job
                void AddJob( const Job & j ) 
                {
                    m_queue.Set(j);
                }
            private:
                JobQueue m_queue ;
                std::vector<std::thread*> m_threads;
        };

    }
}
#endif //__COMMON_MULTITHREAD_MULTITHREAD_H__
