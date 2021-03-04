cmake_minimum_required(VERSION 3.5)

project(fl4 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CXXFLAGS -fopenmp)

#find_package(Boost 1.71.0)

find_package(VTK COMPONENTS
  vtkCommonColor
  vtkCommonCore
  vtkFiltersSources
  vtkFiltersGeometry
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkViewsQt
  vtkIOGeometry
  vtkRenderingAnnotation
  vtkRenderingLabel
  vtkIOImage
  QUIET
)
if (NOT VTK_FOUND)
  message("${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5Charts)
add_executable(fl4
    #    gui/main_test.cpp
    gui/main.cpp
    gui/windows/mainwindow.cpp
    gui/windows/mainwindow.h
    gui/windows/mainwindow.ui
    gui/icons.qrc
#    gui/preferences_gui.h

    #property browser
    gui/propertybrowser/qtvariantproperty.h
    gui/propertybrowser/qtbuttonpropertybrowser.cpp
    gui/propertybrowser/qtgroupboxpropertybrowser.cpp
    gui/propertybrowser/qtgroupboxpropertybrowser.h
    gui/propertybrowser/qteditorfactory.cpp
    gui/propertybrowser/qteditorfactory.h
    gui/propertybrowser/qtpropertybrowserutils.cpp
    gui/propertybrowser/qtvariantproperty.cpp
    gui/propertybrowser/qtpropertymanager.cpp
    gui/propertybrowser/qtpropertymanager.h
    gui/propertybrowser/qtpropertybrowserutils_p.h
    gui/propertybrowser/qtpropertybrowser.qrc
    gui/propertybrowser/qttreepropertybrowser.h
    gui/propertybrowser/qttreepropertybrowser.cpp
    gui/propertybrowser/qtpropertybrowser.h
    gui/propertybrowser/qtbuttonpropertybrowser.h
    gui/propertybrowser/qtpropertybrowser.cpp
    gui/objectpropertybrowser/objectpropertybrowser.cpp
    gui/objectpropertybrowser/objectpropertybrowser.h

    gui/backgroundworker.cpp
    gui/backgroundworker.h

    simulation/modelcontrollerinterface.cpp
    simulation/modelcontrollerinterface.h

    simulation/modelcontroller.cpp
    simulation/modelcontroller.h

    geometry/node.cpp
    geometry/node.h

    geometry/element.cpp
    geometry/element.h

    geometry/mesh.cpp
    geometry/mesh.h

    simulation/parameters_sim.h
    simulation/model.cpp
    simulation/model.h

    simulation/equationofmotionsolver.cpp
    simulation/equationofmotionsolver.h

    gui/interactor/SpecialSelector2D.cpp
    gui/interactor/SpecialSelector2D.h

#[[
    geometry/edge.cpp
    geometry/edge.h

    geometry/floevisualization.h
    geometry/floevisualization.cpp
#]]

#[[
    simulation/linearsystem.cpp
    simulation/linearsystem.h
    simulation/serializer.cpp
    simulation/serializer.h
    simulation/SimpleObjectPool.h

    geometry/solid3d.cpp
    geometry/solid3d.h



    simulation/frameinfo.h
#]]
)


target_include_directories(fl4 PUBLIC gui geometry simulation gui/propertybrowser)
target_include_directories(fl4 PUBLIC gui gui/objectpropertybrowser gui/windows)
target_include_directories(fl4 PUBLIC gui gui/interactor)

# some paths need to be set manually
target_include_directories(fl4 PUBLIC /home/s2/Projects-libraries/eigen)
target_include_directories(fl4 PUBLIC /home/s2/Projects-libraries/oneTBB/include/)
target_include_directories(fl4 PUBLIC /home/s2/Projects-libraries/oneTBB/include/tbb/)
target_include_directories(fl4 PUBLIC /home/s2/Projects-libraries/vtk/Common/Core/)

target_link_directories(fl4 PUBLIC /usr/local/lib)
target_link_directories(fl4 PUBLIC /home/s2/Projects-libraries/oneTBB/build/linux_intel64_gcc_cc9.3.0_libc2.31_kernel5.4.0_release)

target_link_libraries(fl4 PRIVATE Qt5::Widgets Qt5::Charts Qt5::Core)
target_link_libraries(fl4 PRIVATE ${VTK_LIBRARIES})
target_link_libraries(fl4 PUBLIC gmsh)
target_link_libraries(fl4 PUBLIC tbb)
target_link_libraries(fl4 PUBLIC gomp pthread m dl)

#MOSEK
target_include_directories(fl4 PUBLIC "/home/s2/mosek/9.2/tools/platform/linux64x86/h")
target_link_directories(fl4 PUBLIC "/home/s2/mosek/9.2/tools/platform/linux64x86/bin")
target_link_libraries(fl4 PRIVATE mosek64)


