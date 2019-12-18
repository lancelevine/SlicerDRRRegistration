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

// FooBar Widgets includes
#include "qSlicerVertRegFooBarWidget.h"
#include "ui_qSlicerVertRegFooBarWidget.h"

#include "vtkMRMLScene.h"


//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_VertReg
class qSlicerVertRegFooBarWidgetPrivate
  : public Ui_qSlicerVertRegFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerVertRegFooBarWidget);
protected:
  qSlicerVertRegFooBarWidget* const q_ptr;

public:
  qSlicerVertRegFooBarWidgetPrivate(
    qSlicerVertRegFooBarWidget& object);
  virtual void setupUi(qSlicerVertRegFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerVertRegFooBarWidgetPrivate
::qSlicerVertRegFooBarWidgetPrivate(
  qSlicerVertRegFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerVertRegFooBarWidgetPrivate
::setupUi(qSlicerVertRegFooBarWidget* widget)
{
  this->Ui_qSlicerVertRegFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerVertRegFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerVertRegFooBarWidget
::qSlicerVertRegFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerVertRegFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerVertRegFooBarWidget);
  d->setupUi(this);

  //connect(d->InputVolumeSelector, SIGNAL(currentNodeChanged(vtkMRMLNode*)),
	//  this, SLOT(setInputVolume(vtkMRMLNode*)));

  //connect(d->InputVolumeSelector, SIGNAL(mrmlSceneChanged(vtkMRMLScene*)),
	//  SLOT(d->InputVolumeSelector->setMRMLScene(vtkMRMLScene*)));

  //scene = new vtkMRMLScene();
  //scene->
  //d->InputVolumeSelector->setMRMLScene(qSlicer)
}

//-----------------------------------------------------------------------------
qSlicerVertRegFooBarWidget
::~qSlicerVertRegFooBarWidget()
{
}


//-----------------------------------------------------------------------------
void qSlicerVertRegFooBarWidget
::setMRMLScene(vtkMRMLScene* scene)
{
	//this->Superclass::setMRMLScene(scene);
}