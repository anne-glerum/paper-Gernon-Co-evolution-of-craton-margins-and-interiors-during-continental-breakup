#include <aspect/simulator.h>

/*
 * Launch the following function when this plugin is created. Launch ASPECT
 * twice to test checkpoint/resume and then abort the outer ASPECT run.
 */
int f()
{
  std::cout << "* starting from beginning:" << std::endl;

  // call ASPECT with "--" and pipe an existing input file into it.
  int ret;

  ret = system ("cd output-resume_free_surface ; "
                "(cat " ASPECT_SOURCE_DIR "/tests/resume_free_surface.prm "
                " ; "
                " echo 'set Output directory = output1.tmp' "
                " ; "
                " rm -rf output1.tmp ; mkdir output1.tmp "
                ") "
                "| ../../aspect -- >/dev/null ");

  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  ret = system ("cd output-resume_free_surface ; "
                " rm -rf output2.tmp ; cp -r output1.tmp output2.tmp ;"
                " rm -f output1.tmp/log.txt; "
                " cp output1.tmp/restart.mesh.old output1.tmp/restart.mesh;"
                " cp output1.tmp/restart.mesh.info.old output1.tmp/restart.mesh.info;"
                " cp output1.tmp/restart.resume.z.old output1.tmp/restart.resume.z;");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;


  std::cout << "* now resuming:" << std::endl;
  ret = system ("cd output-resume_free_surface ; "
                "(cat " ASPECT_SOURCE_DIR "/tests/resume_free_surface.prm "
                " ; "
                " echo 'set Output directory = output1.tmp' "
                " ; "
                " echo 'set Resume computation = true' "
                ") "
                "| ../../aspect -- >/dev/null");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  std::cout << "* now comparing:" << std::endl;

  ret = system ("cd output-resume_free_surface ; "
                " cp output1.tmp/solution/solution-00004.0000.gnuplot solution/solution-00004.0000.gnuplot1 ;"
                " cp output2.tmp/solution/solution-00004.0000.gnuplot solution/solution-00004.0000.gnuplot2 ;"
                " cp output1.tmp/log.txt log.txt1 ;"
                " cp output2.tmp/log.txt log.txt2 ;"
                " cp output1.tmp/statistics statistics1 ;"
                " cp output2.tmp/statistics statistics2");

  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  ret = system ("cd output-resume_free_surface ; "
                "diff -c output?.tmp/restart.resume.z ; "
                "diff -c output?.tmp/restart.mesh ; "
                "");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  // abort current process:
  exit (0);
  return 42;
}


// run this function by initializing a global variable by it
int i = f();
