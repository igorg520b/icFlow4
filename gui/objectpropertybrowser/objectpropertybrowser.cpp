#include "objectpropertybrowser.h"
#include <QDebug>

ObjectPropertyBrowser::ObjectPropertyBrowser(QWidget* parent) : QtTreePropertyBrowser(parent)
{
    variantManager = new QtVariantPropertyManager(this);
    setFactoryForManager(variantManager, new QtVariantEditorFactory);
}

void ObjectPropertyBrowser::setActiveObject(QObject *obj)
{
    clear();
    variantManager->clear();
    propertyMap.clear();
    if(currentlyConnectedObject) currentlyConnectedObject->disconnect(this);
    currentlyConnectedObject = obj;
    if(!obj) return;

    int count = obj->metaObject()->propertyCount();
    for(int i=1; i< count;i++) {

        QMetaProperty mp = obj->metaObject()->property(i);
        QtVariantProperty *property = variantManager->addProperty(mp.type(),mp.name());
        if(!property) continue;

        property->setEnabled(mp.isWritable());
        propertyMap[property] = mp.name();
        addProperty(property);
    }
    connect(obj, SIGNAL(propertyChanged()), this, SLOT(objectUpdated()));
    objectUpdated();
}

void ObjectPropertyBrowser::valueChanged(QtProperty *property, const QVariant &value)
{
    currentlyConnectedObject->setProperty(propertyMap[property], value);
    objectUpdated();
}

void ObjectPropertyBrowser::objectUpdated()
{
    disconnect(variantManager, SIGNAL(valueChanged(QtProperty*, QVariant)), this, SLOT(valueChanged(QtProperty*, QVariant)));
    QMapIterator<QtProperty*, const char*> i(propertyMap);
    while(i.hasNext()) {
        i.next();
        variantManager->setValue(i.key(), currentlyConnectedObject->property(i.value()));
    }
    connect(variantManager, SIGNAL(valueChanged(QtProperty*, QVariant)), this, SLOT(valueChanged(QtProperty*, QVariant)));
}


