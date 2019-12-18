/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

#ifndef __qSlicerVertRegFooBarWidget_h
#define __qSlicerVertRegFooBarWidget_h

// Qt includes
#include <QWidget>

// FooBar Widgets includes
#include "qSlicerVertRegModuleWidgetsExport.h"

#include "vtkMRMLScene.h"


class qSlicerVertRegFooBarWidgetPrivate;

/// \ingroup Slicer_QtModules_VertReg
class Q_SLICER_MODULE_VERTREG_WIDGETS_EXPORT qSlicerVertRegFooBarWidget
  : public QWidget
{
  Q_OBJECT
public:
  typedef QWidget Superclass;
  qSlicerVertRegFooBarWidget(QWidget *parent=0);
  virtual ~qSlicerVertRegFooBarWidget();

public slots:

  void setMRMLScene(vtkMRMLScene* scene);

protected slots:

protected:
  QScopedPointer<qSlicerVertRegFooBarWidgetPrivate> d_ptr;
  //vtkMRMLScene scene;

private:
  Q_DECLARE_PRIVATE(qSlicerVertRegFooBarWidget);
  Q_DISABLE_COPY(qSlicerVertRegFooBarWidget);
};

#endif
