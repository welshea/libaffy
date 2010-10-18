#ifndef __UTIL_LOG_H__
#define __UTIL_LOG_H__

/*********************************************************************
 * 
 * Logging utilities.
 * The method used is to define generic macros that capture some extra
 * information, etc. and then pass them off to actual implementations.
 * These implementations are part of separate libraries (libtxtlog and
 * libwxlog), so that different types of implementations are provided.
 *********************************************************************/

/* Logging utilities */
#ifdef __cplusplus
extern "C"
{
#endif

  /* Add #ifdef's here eventually to select between logging types */
  #define LIBUTILS_MAX_PB_DEPTH 16
  typedef struct libutils_pb_state_s
  {
    unsigned int depth;
    unsigned int cur_ticks[LIBUTILS_MAX_PB_DEPTH];
    unsigned int tick_interval[LIBUTILS_MAX_PB_DEPTH];
    unsigned int max[LIBUTILS_MAX_PB_DEPTH];
    void *ptr[LIBUTILS_MAX_PB_DEPTH];
  } LIBUTILS_PB_STATE;

  #define LIBUTILS_STR1(x) #x
  #define LIBUTILS_STR(x) LIBUTILS_STR1(x)

  #define _UTILS_F __FILE__
  #define _UTILS_L LIBUTILS_STR(__LINE__)
  #define _UTILS_FL _UTILS_F ":" _UTILS_L ":"

  void die(const char *msg, ...);
  void debug(const char *msg, ...);
  void warn(const char *msg, ...);
  void info(const char *msg, ...);
  void status(const char *msg, ...);

/* Progress bar utilities */
  void pb_init(LIBUTILS_PB_STATE *pbs);
  void pb_cleanup(LIBUTILS_PB_STATE *pbs);
  void pb_begin(LIBUTILS_PB_STATE *pbs, unsigned int max, char *title, ...);
  void pb_tick(LIBUTILS_PB_STATE *pbs, unsigned int tick_sz, char *msg,...);
  void pb_msg(LIBUTILS_PB_STATE *pbs, char *msg, ...);
  void pb_finish(LIBUTILS_PB_STATE *pbs, char *msg, ...);
#ifdef __cplusplus
}
#endif

#endif
