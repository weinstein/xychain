#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <memory>
#include <pthread.h>
#include <queue>
using std::queue;

class Notification;

class ThreadPool {
 public:
   ThreadPool(int nthreads);
   ~ThreadPool();

   void AddJob(void (*func)(void*), void* arg);
   void AddJobWithNotification(void (*func)(void*), void* arg, Notification* notify);
   int GetPendingJobs() const;
   int GetNthreads() const;
   void Wait();

 private:
   struct Job {
      void (*func)(void*);
      void* arg;
      Notification* notify;
   };

   int GetPendingJobsInternal() const;
   static void* ThreadLoop(void* arg);
   void Shutdown();

   mutable pthread_mutex_t fMu;
   mutable pthread_cond_t fCond;
   const int fNthreads;
   int fNpending;
   std::unique_ptr<pthread_t[]> fTIDs;

   queue<Job> fJobQueue;

   bool fIsShutdown;
};

#endif   
