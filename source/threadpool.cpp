#include <memory>
#include <pthread.h>
#include <queue>

#include "threadpool.h"
#include "mutex.h"
#include "notification.h"

ThreadPool::ThreadPool(int nthreads) : fNthreads(nthreads),
                                       fNpending(0),
                                       fIsShutdown(false),
                                       fTIDs(new pthread_t[fNthreads]) {
   pthread_mutex_init(&fMu, NULL);
   pthread_cond_init(&fCond, NULL);

   for (int i = 0; i < fNthreads; i++) {
      pthread_create(&fTIDs[i], NULL, ThreadLoop, (void*)this);
   }
}

ThreadPool::~ThreadPool() {
   Shutdown();
}

void ThreadPool::AddJobWithNotification(void (*func)(void*), void* arg, Notification* notify) {
   ThreadPool::Job job = {func, arg, notify};

   {
      MutexLock ml(&fMu);
      fJobQueue.push(job);
   }

   pthread_cond_signal(&fCond);
}

void ThreadPool::AddJob(void (*func)(void*), void* arg) {
   AddJobWithNotification(func, arg, NULL);
}

int ThreadPool::GetPendingJobs() const {
   MutexLock ml(&fMu);
   return GetPendingJobsInternal();
}

int ThreadPool::GetPendingJobsInternal() const {
   return fJobQueue.size() + fNpending;
}

int ThreadPool::GetNthreads() const {
   return fNthreads;
}

void ThreadPool::Wait() {
   MutexLock ml(&fMu);
   while (GetPendingJobsInternal() > 0)
      ml.CondWait(&fCond);
}

void ThreadPool::Shutdown() {
   {
      MutexLock ml(&fMu);
      fIsShutdown = true;
      while (!fJobQueue.empty())
         fJobQueue.pop();
   }
   pthread_cond_broadcast(&fCond);

   void* status;
   for (int i = 0; i < fNthreads; i++) {
      pthread_join(fTIDs[i], &status);
   }

   pthread_mutex_destroy(&fMu);
   pthread_cond_destroy(&fCond);
}

void* ThreadPool::ThreadLoop(void* arg) {
   ThreadPool* self = (ThreadPool*)arg;
   while (1) {
      ThreadPool::Job job;

      {
         MutexLock ml(&self->fMu);
         while (self->fJobQueue.empty() && !self->fIsShutdown) {
            ml.CondWait(&self->fCond);
         }
         if (self->fIsShutdown) {
            break;
         }

         job = self->fJobQueue.front();
         self->fJobQueue.pop();
         ++self->fNpending;
      }

      (job.func)(job.arg);
      if (job.notify) {
         job.notify->Notify();
      }

      {
         MutexLock ml(&self->fMu);
         --self->fNpending;
      }
      pthread_cond_broadcast(&self->fCond);
   }

   pthread_exit(0);
}
