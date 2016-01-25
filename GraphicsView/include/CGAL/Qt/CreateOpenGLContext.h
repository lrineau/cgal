// Copyright (c) 2015  GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent Rineau and Maxime Gimeno
#ifndef CGAL_QT_CREATE_OPENGL_CONTEXT_H
#define CGAL_QT_CREATE_OPENGL_CONTEXT_H

#include <CGAL/assertions.h>
#include <QOpenGLContext>
#include <QGLContext>
#include <QString>
#include <QSurfaceFormat>

namespace CGAL{
namespace Qt{
inline QGLContext* createOpenGLContext()
{
    QOpenGLContext *context = new QOpenGLContext();
    QSurfaceFormat format;
    format.setVersion(2,1);
    format.setProfile(QSurfaceFormat::CompatibilityProfile);
    context->setFormat(format);
    context->create();
    QGLContext* old_context = QGLContext::fromOpenGLContext(context);
    if(!old_context->isValid()) {
      CGAL_error_msg("Cannot create the OpenGL context!");
    } else
    {
      QSurfaceFormat format = old_context->contextHandle()->format();
      QGLFormat old_format = old_context->format();
      std::cerr << "QGLcontext created, "
                << "context: " << (void*)context << "\n"
                << "old_context: " << (void*)old_context << "\n";
      std::cerr << "old_context->contextHandle(): "
                << (void*)old_context->contextHandle() << "\n";
      std::cerr << "old_context->isValid(): " << old_context->isValid() << "\n";
      std::cerr << qPrintable(QString("  format: OpenGL %1.%2%3\n")
                              .arg(format.majorVersion())
                              .arg(format.minorVersion())
                              .arg(format.profile() == QSurfaceFormat::CoreProfile ?
                                   QString(" core profile") :
                                   format.profile() == QSurfaceFormat::CompatibilityProfile ?
                                   QString(" compatibility profile") :
                                   QString("")));
      std::cerr << qPrintable(QString("  old format: OpenGL %1.%2%3\n")
                              .arg(old_format.majorVersion())
                              .arg(old_format.minorVersion())
                              .arg(old_format.profile() == QGLFormat::CoreProfile ?
                                   QString(" core profile") :
                                   old_format.profile() == QGLFormat::CompatibilityProfile ?
                                   QString(" compatibility profile") :
                                   QString("")));
    }
    return old_context;
}
} // namespace Qt
} // namespace CGAL
#endif
