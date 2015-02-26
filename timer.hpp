/* 
 * File:   timer.h
 * Author: aksel
 *
 * Created on 12. Oktober 2014, 15:11
 */

#ifndef TIMER_H
#define	TIMER_H

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#define POSIX
#endif


#ifdef POSIX

#include <sys/time.h>
#include <ctime>

namespace nanolens
{
  namespace util
  {
    /// A timer to stop time with nanoseconds precision. Only available on
    /// POSIX systems.

    class timer
    {
      timespec start_;
      timespec stop_;
      bool is_running_;
    public:
      /// Construct object

      timer()
      : is_running_(false)
      {
      }

      /// start the timer

      void start()
      {
        clock_gettime(CLOCK_MONOTONIC, &start_);
        is_running_ = true;
      }

      /// @return whether the timer is currently running

      bool is_running() const
      {
        return is_running_;
      }

      /// Stops the timer.
      /// @return The duration since the call to \c timer::start() in seconds
      /// with nanoseconds of precision

      double stop()
      {
        if (!is_running_)
          return 0.0;

        clock_gettime(CLOCK_MONOTONIC, &stop_);
        double t = stop_.tv_sec - start_.tv_sec;
        t += static_cast<double> (stop_.tv_nsec - start_.tv_nsec) * 1e-9;

        is_running_ = false;
        return t;
      }
    };

    /// A timer for MPI-parallel programs. The \c start() and \c stop() methods
    /// wait until all processes have reached the same position in the code
    /// before starting or stopping the timer.

    class mpi_synced_timer
    {
      MPI_Comm comm_;
      timer t_;
    public:
      /// Initializes the timer
      /// @param communicator The MPI communicator that is used by all processes
      /// that shall be timed. If \c communicator is \c MPI_COMM_NULL, the
      /// object behaves like a normal, non-parallel \c timer object.

      mpi_synced_timer(MPI_Comm communicator)
      : comm_(communicator)
      {
      }

      /// @return whether the timer is currently running

      bool is_running() const
      {
        return t_.is_running();
      }

      /// start the timer. Collective on the supplied communicator, unless
      /// it was MPI_COMM_NULL.

      void start()
      {
        if (comm_ != MPI_COMM_NULL)
          MPI_Barrier(comm_);
        t_.start();
      }

      /// Stops the timer. Collective on the supplied communicator, unless
      /// it was MPI_COMM_NULL.
      /// @return The duration since the call to \c timer::start() in seconds
      /// with nanoseconds of precision

      double stop()
      {
        if (comm_ != MPI_COMM_NULL)
          MPI_Barrier(comm_);
        return t_.stop();
      }
    };

  }
}

#endif //POSIX

#endif	/* TIMER_H */

