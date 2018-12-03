#include "Polyhedron_demo.h"
#include "MainWindow.h"
#include <QMessageBox>
#include <CGAL/Qt/resources.h>
#include <stdexcept>

#include <QCommandLineParser>
#include <QCommandLineOption>
#include <QSurfaceFormat>
#include <QOpenGLContext>

bool use_value_outside;
int interpolation_method;
int interpolation_width;
double gaussian_epsilon;

struct Polyhedron_demo_impl {
  bool catch_exceptions;
  QScopedPointer<MainWindow> mainWindow;

  Polyhedron_demo_impl() : catch_exceptions(true) {}
}; // end struct Polyhedron_demo_impl

Polyhedron_demo::Polyhedron_demo(int& argc, char **argv,
                                 QString application_name,
                                 QString main_window_title)
  : QApplication(argc, argv)
  , d_ptr_is_initialized(false)
  , d_ptr(new Polyhedron_demo_impl)
{
  d_ptr_is_initialized = true;
  std::cerr.precision(17);
  std::cout.precision(17);
  std::clog.precision(17);

  //for windows
#if (QT_VERSION >= QT_VERSION_CHECK(5, 3, 0))
  this->setAttribute(Qt::AA_UseDesktopOpenGL);
#endif

  // Import resources from libCGAL (Qt5).
  CGAL_QT_INIT_RESOURCES;

  this->setOrganizationDomain("geometryfactory.com");
  this->setOrganizationName("GeometryFactory");
  this->setApplicationName(application_name);

  QCommandLineParser parser;
  parser.addHelpOption();

  QCommandLineOption use_meta("use-meta",
                              tr("Use the [Meta] key to move frames, instead of [Tab]."));
  parser.addOption(use_meta);
  QCommandLineOption no_try_catch("no-try-catch",
                                  tr("Do not catch uncaught exceptions."));
  parser.addOption(no_try_catch);
  QCommandLineOption debug_scripts("debug-scripts",
                                   tr("Use the scripts debugger."));
  parser.addOption(debug_scripts);
  QCommandLineOption no_debug_scripts("no-debug-scripts",
                                   tr("Do not use the scripts debugger."));
  parser.addOption(no_debug_scripts);
  QCommandLineOption no_autostart("no-autostart",
                                  tr("Ignore the autostart.js file, if any."));
  parser.addOption(no_autostart);
  QCommandLineOption verbose("verbose",
                                   tr("Print the paths explored byt the application searching for plugins."));
  parser.addOption(verbose);
  QCommandLineOption do_not_use_value_outside("do_not_use_value_outside",
                                              tr("Mesh_3: do not use 'value_outside'"));
  parser.addOption(do_not_use_value_outside);
  QCommandLineOption interpolation_method
    ("interpolation_method",
     tr("Mesh_3: interpolation method (linear: 0, cubic: 1, gaussian: 2)"),
     tr("<method>"), "0");
  parser.addOption(interpolation_method);
  QCommandLineOption interpolation_width("interpolation_width",
                                         tr("Mesh_3: interpolation width"),
                                         tr("<width>"), "2");
  parser.addOption(interpolation_width);
  QCommandLineOption gaussian_epsilon("gaussian_epsilon",
                                      tr("Mesh_3: Gaussian epsilon"),
                                      tr("<epsilon>"), "-1.");
  parser.addOption(gaussian_epsilon);
  QCommandLineOption old("old",
    tr("Force OpenGL 2.1 context."));
  parser.addOption(old);
  parser.addPositionalArgument("files", tr("Files to open"), "[files...]");
  parser.process(*this);
  d_ptr->mainWindow.reset(new MainWindow(parser.isSet(verbose)));
  MainWindow& mainWindow = *d_ptr->mainWindow;
  
  ::use_value_outside = ! parser.isSet(do_not_use_value_outside);
  ::interpolation_method = parser.value(interpolation_method).toInt();
  ::interpolation_width = parser.value(interpolation_width).toInt();
  ::gaussian_epsilon = parser.value(gaussian_epsilon).toDouble();
  mainWindow.setWindowTitle(main_window_title);
  mainWindow.show();

  // On Apple, the first time the application is launched, the menus are unclicable, and
  // the only way you can fix it is to unfocus and re-focus the application.
  // This is a hack that makes the application lose the focus after it is started, to force the user
  // to re-focus it. (source : http://www.alecjacobson.com/weblog/?p=3910)
#ifdef __APPLE__
    system("osascript -e 'tell application \"System Events\" "
      "to keystroke tab using {command down, shift down}'");
#endif
  if(parser.isSet(use_meta)) {
    mainWindow.setAddKeyFrameKeyboardModifiers(::Qt::MetaModifier);
  }
  if(parser.isSet(no_try_catch)) {
    this->do_not_catch_exceptions();
  }
#ifdef QT_SCRIPTTOOLS_LIB
  if(parser.isSet(debug_scripts)) {
    mainWindow.enableScriptDebugger();
  }
  if(parser.isSet(no_debug_scripts)) {
    mainWindow.enableScriptDebugger(false);
  }
#else
  if(parser.isSet(debug_scripts)) {
    std::cerr << "Qt Script Tools have not been configured!";
    abort();
  }
#endif

  mainWindow.loadScript(":/cgal/Polyhedron_3/javascript/lib.js");
  QFileInfo autostart_js("autostart.js");
  if(!parser.isSet(no_autostart) && autostart_js.exists()) {
    mainWindow.loadScript(autostart_js);
  }
  Q_FOREACH(QString filename, parser.positionalArguments()) {
    mainWindow.open(filename);
  }

}

Polyhedron_demo::~Polyhedron_demo() {}

void Polyhedron_demo::do_not_catch_exceptions() {
  d_ptr->catch_exceptions = false;
  setProperty("no-try-catch", true);
}

bool Polyhedron_demo::notify(QObject* receiver, QEvent* event)
{
  if(!d_ptr_is_initialized || !d_ptr->catch_exceptions)
    return QApplication::notify(receiver, event);
  else try {
      return QApplication::notify(receiver, event);
    } catch (std::exception &e) {
      // find the mainwindow to spawn an error message
      Q_FOREACH (QWidget *widget, QApplication::topLevelWidgets()) {
        if(MainWindow* mw = qobject_cast<MainWindow*>(widget)) {
          QMessageBox::critical(mw,
                                tr("Unhandled exception"),
                                tr("<p>Unhandled exception:</p>\n"
                                   "<pre>%1</pre>").arg(e.what()));
          break;
        }
      }
      QApplication::restoreOverrideCursor();
    } catch (...) {
      qFatal("Unknown exception encountered. Aborting.");
    }
  return false;
}

int Polyhedron_demo::try_exec()
{
  // A Qt Script may have closed the main window.
  // The following loop launch app.exec() only if the main window is visible.
  if(d_ptr->mainWindow->isVisible()) {
    return this->exec();
  }
  return 0;
}
